"""A small YAML reader for celSetup input files.

It reads the subset celSetup files need: nested mappings by indentation, scalars (numbers, yes/no,
true/false, strings), inline lists [a, b] and comments. Always this reader,
never PyYAML (which reads 1e-3 as a string), so the tools behave the same
everywhere and run on a bare OpenFOAM image with only numpy.
"""


def _scalar(s):
    s = s.strip()
    if not s:
        return None
    if s[0] in "\"'" and s[-1] == s[0]:
        return s[1:-1]
    if s.startswith("[") and s.endswith("]"):
        inner = s[1:-1].strip()
        return [_scalar(x) for x in inner.split(",")] if inner else []
    low = s.lower()
    if low in ("true", "yes", "on"):
        return True
    if low in ("false", "no", "off"):
        return False
    if low in ("null", "~"):
        return None
    try:
        return int(s)
    except ValueError:
        pass
    try:
        return float(s)
    except ValueError:
        return s


def _strip_comment(line):
    out, q = [], None
    for ch in line:
        if q:
            if ch == q:
                q = None
        elif ch in "\"'":
            q = ch
        elif ch == "#":
            break
        out.append(ch)
    return "".join(out).rstrip()


def loads(text):
    root = {}
    stack = [(-1, root)]
    for n, raw in enumerate(text.splitlines(), 1):
        line = _strip_comment(raw)
        if not line.strip():
            continue
        ind = len(line) - len(line.lstrip(" "))
        if "\t" in line[:ind + 1]:
            raise ValueError(f"line {n}: use spaces, not tabs, for indentation")
        key, sep, val = line.strip().partition(":")
        if not sep:
            raise ValueError(f"line {n}: expected 'key: value', got {raw.strip()!r}")
        while stack[-1][0] >= ind:
            stack.pop()
        parent = stack[-1][1]
        key = key.strip()
        if val.strip():
            parent[key] = _scalar(val)
        else:
            parent[key] = {}
            stack.append((ind, parent[key]))
    return root


def load(path):
    with open(path) as f:
        return loads(f.read())

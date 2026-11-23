/* celFoam user manual: builds the side menu, the on-page contents and the
   previous/next links from one list (NAV), and renders the equations with the
   KaTeX copy in assets/katex. Pages stay readable without JavaScript; this only
   adds navigation and typesetting. tests/test_manual.sh checks that NAV and the
   files in docs/manual agree. */

var NAV = [
    { group: "Start here", pages: [
        ["index.html", "Home"],
        ["install.html", "Installation"],
        ["quickstart.html", "Quick start"]
    ]},
    { group: "Using celFoam", pages: [
        ["case-setup.html", "Case setup reference"],
        ["meshing.html", "Meshing guide"]
    ]},
    { group: "Theory", pages: [
        ["theory-equations.html", "Equations and numerics"],
        ["theory-materials.html", "Material models"],
        ["theory-bodies.html", "Bodies and contact"]
    ]},
    { group: "Validation", pages: [
        ["validation.html", "Validation cases"]
    ]},
    { group: "Reference", pages: [
        ["faq.html", "FAQ and known issues"],
        ["references.html", "References"]
    ]}
];

(function ()
{
    "use strict";

    var here = location.pathname.split("/").pop() || "index.html";

    function el(tag, attrs, text)
    {
        var e = document.createElement(tag);
        for (var k in attrs) { e.setAttribute(k, attrs[k]); }
        if (text) { e.textContent = text; }
        return e;
    }

    // ---- theme: follow the system, with a manual override kept per browser ----
    function storedTheme()
    {
        try { return localStorage.getItem("celManualTheme"); } catch (e) { return null; }
    }

    function applyTheme(t)
    {
        if (t) { document.documentElement.setAttribute("data-theme", t); }
        else { document.documentElement.removeAttribute("data-theme"); }
    }

    applyTheme(storedTheme());

    function build()
    {
        var flat = [];
        NAV.forEach(function (g) { g.pages.forEach(function (p) { flat.push(p); }); });

        // ---- side menu ----
        var side = el("nav", { "class": "side", "aria-label": "Manual contents" });
        var brand = el("a", { "class": "brand", href: "index.html" }, "celFoam manual");
        side.appendChild(brand);
        side.appendChild(el("p", { "class": "tagline" },
            "Explicit Coupled Eulerian-Lagrangian solver for bulk solids, on OpenFOAM v2512"));

        NAV.forEach(function (g)
        {
            side.appendChild(el("div", { "class": "group" }, g.group));
            var ul = el("ul");
            g.pages.forEach(function (p)
            {
                var li = el("li");
                var a = el("a", { "class": "page" + (p[0] === here ? " current" : ""), href: p[0] }, p[1]);
                li.appendChild(a);
                if (p[0] === here)
                {
                    var heads = document.querySelectorAll("main h2[id]");
                    if (heads.length)
                    {
                        var toc = el("div", { "class": "toc" });
                        heads.forEach(function (h)
                        {
                            toc.appendChild(el("a", { href: "#" + h.id }, h.textContent));
                        });
                        li.appendChild(toc);
                    }
                }
                ul.appendChild(li);
            });
            side.appendChild(ul);
        });

        var tbtn = el("button", { "class": "theme-toggle", type: "button" }, "Light / dark");
        tbtn.addEventListener("click", function ()
        {
            var dark = document.documentElement.getAttribute("data-theme") === "dark"
                || (!document.documentElement.getAttribute("data-theme")
                    && window.matchMedia("(prefers-color-scheme: dark)").matches);
            var next = dark ? "light" : "dark";
            applyTheme(next);
            try { localStorage.setItem("celManualTheme", next); } catch (e) { /* private window */ }
        });
        side.appendChild(tbtn);

        // ---- top bar for narrow screens ----
        var bar = el("div", { "class": "topbar" });
        var mbtn = el("button", { type: "button", "aria-label": "Open the menu" }, "Menu");
        mbtn.addEventListener("click", function () { document.body.classList.toggle("menu-open"); });
        bar.appendChild(mbtn);
        bar.appendChild(el("a", { href: "index.html" }, "celFoam manual"));

        document.body.insertBefore(side, document.body.firstChild);
        document.body.insertBefore(bar, document.body.firstChild);

        side.addEventListener("click", function (ev)
        {
            if (ev.target.tagName === "A") { document.body.classList.remove("menu-open"); }
        });

        // ---- previous / next ----
        var i = -1;
        flat.forEach(function (p, k) { if (p[0] === here) { i = k; } });
        var main = document.querySelector("main");
        if (main && i >= 0)
        {
            var pager = el("div", { "class": "pager" });
            if (i > 0)
            {
                var pa = el("a", { "class": "prev", href: flat[i - 1][0] });
                pa.appendChild(el("small", {}, "Previous"));
                pa.appendChild(document.createTextNode(flat[i - 1][1]));
                pager.appendChild(pa);
            }
            if (i < flat.length - 1)
            {
                var na = el("a", { "class": "next", href: flat[i + 1][0] });
                na.appendChild(el("small", {}, "Next"));
                na.appendChild(document.createTextNode(flat[i + 1][1]));
                pager.appendChild(na);
            }
            main.appendChild(pager);
        }

        // ---- wide tables scroll sideways on phones ----
        document.querySelectorAll("main table").forEach(function (t)
        {
            if (t.parentNode.classList && t.parentNode.classList.contains("table-wrap")) { return; }
            var w = el("div", { "class": "table-wrap" });
            t.parentNode.insertBefore(w, t);
            w.appendChild(t);
        });

        // ---- equations ----
        if (window.renderMathInElement)
        {
            window.renderMathInElement(document.querySelector("main"), {
                delimiters: [
                    { left: "$$", right: "$$", display: true },
                    { left: "\\(", right: "\\)", display: false }
                ],
                throwOnError: false
            });
        }
    }

    if (document.readyState === "loading") { document.addEventListener("DOMContentLoaded", build); }
    else { build(); }
})();

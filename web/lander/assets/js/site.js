/* popkern site: the only script on the page. No dependencies, no build. */
(function () {
  "use strict";

  /* Copy-to-clipboard on code blocks. */
  document.querySelectorAll(".copy-btn").forEach(function (btn) {
    btn.addEventListener("click", function () {
      var wrap = btn.closest(".code-wrap");
      var code = wrap && wrap.querySelector("code");
      if (!code) return;

      var done = function (ok) {
        btn.textContent = ok ? "copied" : "copy failed";
        btn.dataset.copied = ok ? "1" : "";
        setTimeout(function () {
          btn.textContent = "copy";
          btn.dataset.copied = "";
        }, 1800);
      };

      var text = code.innerText;
      if (navigator.clipboard && window.isSecureContext) {
        navigator.clipboard.writeText(text).then(function () { done(true); }, function () { done(false); });
      } else {
        /* http:// preview fallback. The clipboard API is secure-context only. */
        var ta = document.createElement("textarea");
        ta.value = text;
        ta.style.cssText = "position:fixed;opacity:0";
        document.body.appendChild(ta);
        ta.select();
        var ok = false;
        try { ok = document.execCommand("copy"); } catch (e) { ok = false; }
        document.body.removeChild(ta);
        done(ok);
      }
    });
  });

  /* Footer year. */
  document.querySelectorAll("[data-year]").forEach(function (el) {
    el.textContent = String(new Date().getFullYear());
  });
})();

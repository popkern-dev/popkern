// The remembered colour theme, on <html> before the page paints. A
// classic script in <head> (modules run after the first paint, which
// flashed the default), so it cannot import local-store.js: it is the
// one other reader of its THEME_KEY, read-only and guarded. theme.js
// owns the choice; storage may be blocked (a private window): then the
// default (CAD) shows.
try {
  const theme = localStorage.getItem("popkern.theme");
  if (theme === "cad" || theme === "midnight" || theme === "light") document.documentElement.dataset.theme = theme;
} catch {
  // storage blocked: the default theme
}

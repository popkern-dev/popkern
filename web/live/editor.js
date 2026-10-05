// A dependency-free script editor: a textarea (native keyboard,
// selection, undo, screen readers) with a highlighted mirror behind
// it and a gutter of line numbers. Highlighting is the popkern Rhai
// dialect: keywords, surface verbs, strings, numbers, comments.

const KEYWORDS = /\b(let|const|fn|if|else|for|in|while|loop|return|break|continue|true|false|try|catch|throw|import|export|switch|do|until)\b/g;
const VERBS = /\b(model|open|profile|sketch|extrude|revolve|fillet|chamfer|shell|union|subtract|intersect|translate|rotate|edges|faces|named|xy|xz|yz|xy_at|xz_at|yz_at|on_face|blind|through_all|up_to_surface|mid_plane|full_circle|rect|circle|polygon|path|line_to|arc_to|close|tool|param|param_expr|set_param|feature|describe|assert_watertight|label_edge|label_face|print|p|expr|conic|linear|by_stations|equal|two_distance|distance_angle|axis_at|axis|insert_before|insert_after|reorder|suppress|unsuppress|to_rhai)\b(?=\s*\()/g;

function escape(text) {
  return text.replace(/&/g, "&amp;").replace(/</g, "&lt;").replace(/>/g, "&gt;");
}

/** Highlight one line of source as HTML. */
export function highlight(line) {
  const tokens = [];
  const re = /(\/\/.*$)|("(?:[^"\\]|\\.)*")|(\b\d+(?:\.\d+)?(?:e[+-]?\d+)?\b)|([A-Za-z_][A-Za-z0-9_]*)|(\s+|.)/g;
  let m;
  while ((m = re.exec(line))) {
    if (m[1]) tokens.push(`<span class="tk-comment">${escape(m[1])}</span>`);
    else if (m[2]) tokens.push(`<span class="tk-string">${escape(m[2])}</span>`);
    else if (m[3]) tokens.push(`<span class="tk-number">${escape(m[3])}</span>`);
    else if (m[4]) {
      const word = m[4];
      KEYWORDS.lastIndex = 0;
      VERBS.lastIndex = 0;
      const rest = line.slice(m.index + word.length);
      if (KEYWORDS.test(word)) tokens.push(`<span class="tk-keyword">${word}</span>`);
      else if (/^\s*\(/.test(rest) && VERBS.test(`${word}(`)) tokens.push(`<span class="tk-verb">${word}</span>`);
      else tokens.push(escape(word));
    } else tokens.push(escape(m[5]));
  }
  return tokens.join("") || " ";
}

export class Editor {
  constructor(root, { onChange, onRun }) {
    this.root = root;
    this.onChange = onChange;
    this.onRun = onRun;
    root.classList.add("editor");
    this.gutter = document.createElement("div");
    this.gutter.className = "editor-gutter";
    this.mirror = document.createElement("pre");
    this.mirror.className = "editor-mirror";
    this.mirror.setAttribute("aria-hidden", "true");
    this.textarea = document.createElement("textarea");
    this.textarea.className = "editor-input";
    this.textarea.spellcheck = false;
    this.textarea.setAttribute("aria-label", "Model script");
    this.textarea.setAttribute("autocapitalize", "off");
    this.textarea.setAttribute("autocomplete", "off");
    root.append(this.gutter, this.mirror, this.textarea);
    this.markedLine = null;
    this.textarea.addEventListener("input", () => {
      this.refresh();
      this.onChange(this.value);
    });
    this.textarea.addEventListener("scroll", () => this.sync());
    this.textarea.addEventListener("keydown", (e) => {
      if ((e.ctrlKey || e.metaKey) && e.key === "Enter") {
        e.preventDefault();
        this.onRun();
      } else if (e.key === "Tab") {
        e.preventDefault();
        const { selectionStart: s, selectionEnd: t } = this.textarea;
        this.textarea.setRangeText("    ", s, t, "end");
        this.refresh();
        this.onChange(this.value);
      }
    });
  }

  get value() {
    return this.textarea.value;
  }

  set value(text) {
    this.textarea.value = text;
    this.refresh();
  }

  /** Whether typing is refused (a read-only example). */
  get readOnly() {
    return this.textarea.readOnly;
  }

  set readOnly(value) {
    this.textarea.readOnly = value;
    this.root.classList.toggle("editor-readonly", value);
  }

  /** Mark the history scrubber's step line (1-based); `null` clears.
   *  Scrolls it into view without moving the caret. */
  markStep(line) {
    this.stepLine = line;
    this.refresh();
    if (line != null) {
      const height = this.textarea.scrollHeight / Math.max(1, this.value.split("\n").length);
      this.textarea.scrollTop = Math.max(0, (line - 4) * height);
      this.sync();
    }
  }

  /** Mark a 1-based line (an error position); `null` clears. */
  mark(line) {
    this.markedLine = line;
    this.refresh();
  }

  /** Move the caret to a 1-based line and column. */
  goTo(line, column) {
    const lines = this.value.split("\n");
    let offset = 0;
    for (let i = 0; i < Math.min(line - 1, lines.length); i++) offset += lines[i].length + 1;
    offset += Math.max(0, (column ?? 1) - 1);
    this.textarea.focus();
    this.textarea.setSelectionRange(offset, offset);
  }

  refresh() {
    const lines = this.value.split("\n");
    this.mirror.innerHTML = lines
      .map((l, i) => `<div class="editor-line${this.markedLine === i + 1 ? " editor-line-error" : ""}${this.stepLine === i + 1 ? " editor-line-step" : ""}">${highlight(l)}</div>`)
      .join("");
    this.gutter.innerHTML = lines.map((_, i) => `<div>${i + 1}</div>`).join("");
    this.sync();
  }

  sync() {
    this.mirror.scrollTop = this.textarea.scrollTop;
    this.mirror.scrollLeft = this.textarea.scrollLeft;
    this.gutter.scrollTop = this.textarea.scrollTop;
  }
}

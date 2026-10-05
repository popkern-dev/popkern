// A project folder as a .zip and back (live-projects-and-history.md
// §5.7). Dependency-free: entries are STOREd (no compression; project
// files are small text), with CRC-32 and the standard local headers
// and central directory, so `unzip`, Python's zipfile and every OS
// archive tool read them. Reading accepts STOREd entries only, which is
// what this writer produces; a compressed archive is refused with a
// message rather than misread.

const CRC_TABLE = (() => {
  const table = new Uint32Array(256);
  for (let n = 0; n < 256; n++) {
    let c = n;
    for (let k = 0; k < 8; k++) c = c & 1 ? 0xedb88320 ^ (c >>> 1) : c >>> 1;
    table[n] = c >>> 0;
  }
  return table;
})();

/** CRC-32 (IEEE) of `bytes`. */
export function crc32(bytes) {
  let crc = 0xffffffff;
  for (const b of bytes) crc = CRC_TABLE[(crc ^ b) & 0xff] ^ (crc >>> 8);
  return (crc ^ 0xffffffff) >>> 0;
}

const encoder = new TextEncoder();

/** `files` (path -> Uint8Array or string) as zip bytes. Entries are
 *  written in path order so the same project zips identically. */
export function zipFiles(files) {
  const entries = Object.keys(files)
    .sort()
    .map((path) => {
      const data = typeof files[path] === "string" ? encoder.encode(files[path]) : files[path];
      return { name: encoder.encode(path), data, crc: crc32(data) };
    });
  const parts = [];
  const central = [];
  let offset = 0;
  for (const e of entries) {
    const local = new DataView(new ArrayBuffer(30));
    local.setUint32(0, 0x04034b50, true);
    local.setUint16(4, 20, true); // version needed
    local.setUint16(6, 0x0800, true); // UTF-8 names
    local.setUint16(8, 0, true); // stored
    local.setUint32(14, e.crc, true);
    local.setUint32(18, e.data.length, true);
    local.setUint32(22, e.data.length, true);
    local.setUint16(26, e.name.length, true);
    parts.push(new Uint8Array(local.buffer), e.name, e.data);
    const header = new DataView(new ArrayBuffer(46));
    header.setUint32(0, 0x02014b50, true);
    header.setUint16(4, 20, true);
    header.setUint16(6, 20, true);
    header.setUint16(8, 0x0800, true);
    header.setUint32(16, e.crc, true);
    header.setUint32(20, e.data.length, true);
    header.setUint32(24, e.data.length, true);
    header.setUint16(28, e.name.length, true);
    header.setUint32(42, offset, true);
    central.push(new Uint8Array(header.buffer), e.name);
    offset += 30 + e.name.length + e.data.length;
  }
  const centralSize = central.reduce((n, p) => n + p.length, 0);
  const end = new DataView(new ArrayBuffer(22));
  end.setUint32(0, 0x06054b50, true);
  end.setUint16(8, entries.length, true);
  end.setUint16(10, entries.length, true);
  end.setUint32(12, centralSize, true);
  end.setUint32(16, offset, true);
  const all = [...parts, ...central, new Uint8Array(end.buffer)];
  const out = new Uint8Array(all.reduce((n, p) => n + p.length, 0));
  let at = 0;
  for (const p of all) {
    out.set(p, at);
    at += p.length;
  }
  return out;
}

const decoder = new TextDecoder();

/** The files of zip `bytes` (path -> Uint8Array). Folders are skipped;
 *  a single top folder shared by every entry (an OS "compress folder")
 *  is stripped. Throws an Error naming what it cannot read. */
export function unzipFiles(bytes) {
  const view = new DataView(bytes.buffer, bytes.byteOffset, bytes.byteLength);
  let end = -1;
  for (let i = bytes.length - 22; i >= Math.max(0, bytes.length - 22 - 65535); i--) {
    if (view.getUint32(i, true) === 0x06054b50) {
      end = i;
      break;
    }
  }
  if (end < 0) throw new Error("not a zip archive (no end of central directory)");
  const count = view.getUint16(end + 10, true);
  let at = view.getUint32(end + 16, true);
  const files = {};
  for (let n = 0; n < count; n++) {
    if (view.getUint32(at, true) !== 0x02014b50) throw new Error("damaged zip central directory");
    const method = view.getUint16(at + 10, true);
    const crc = view.getUint32(at + 16, true);
    const size = view.getUint32(at + 20, true);
    const nameLength = view.getUint16(at + 28, true);
    const extraLength = view.getUint16(at + 30, true);
    const commentLength = view.getUint16(at + 32, true);
    const localAt = view.getUint32(at + 42, true);
    const name = decoder.decode(bytes.subarray(at + 46, at + 46 + nameLength));
    at += 46 + nameLength + extraLength + commentLength;
    if (name.endsWith("/")) continue;
    if (method !== 0) throw new Error(`${name} is compressed; export projects from popkern live, or zip them without compression (zip -0)`);
    const dataAt = localAt + 30 + view.getUint16(localAt + 26, true) + view.getUint16(localAt + 28, true);
    const data = bytes.slice(dataAt, dataAt + size);
    if (crc32(data) !== crc) throw new Error(`${name} is damaged (checksum mismatch)`);
    files[name] = data;
  }
  // An archive of the folder itself ("my-bracket/popkern.toml", …):
  // strip the shared top folder.
  const paths = Object.keys(files);
  const top = paths[0]?.includes("/") ? `${paths[0].split("/")[0]}/` : null;
  if (!("popkern.toml" in files) && top && paths.every((p) => p.startsWith(top))) {
    const stripped = {};
    for (const p of paths) stripped[p.slice(top.length)] = files[p];
    return stripped;
  }
  return files;
}

# Vendored browser dependencies

The page loads no third-party code over the network. Each vendored
file is pinned by exact version and content hash; the popkern-live page test
suite (`crates/popkern-live/tests/page.rs`) fails when a hash
changes without this file being updated.

| File | Project | Version | License | SHA-256 |
|---|---|---|---|---|
| `three.module.min.js` | three.js (https://threejs.org) | 0.170.0 (r170) | MIT — Copyright © 2010-2024 three.js authors | `08fd7545d13d2c7fb65ab691530a802dafefd638596501854f267d0fb13c39e7` |
| `mp4-muxer.js` | mp4-muxer (https://github.com/Vanilagy/mp4-muxer) | 5.2.2 | MIT — Copyright (c) 2023 Vanilagy | `d2c4c782f180c86ed30b1f5d9487a34a0d370bf9b4535285734ea187d38f9bb5` |
| `webm-muxer.js` | webm-muxer (https://github.com/Vanilagy/webm-muxer) | 5.1.4 | MIT — Copyright (c) 2022 Vanilagy | `53710543cecc4db7b6cd23fb0d9b308aa3374697f43a6149a8497f0a8547f0cb` |
| `gifenc.esm.js` | gifenc (https://github.com/mattdesl/gifenc) | 1.0.3 | MIT — Copyright (c) 2017 Matt DesLauriers | `3767ab0eef10bb10fd4bacf8c708fa307fa59062541d24a307e4c842b66bba80` |

Source: `https://cdnjs.cloudflare.com/ajax/libs/three.js/0.170.0/three.module.min.js`,
fetched 2026-09-02. The MIT license text:

```
Copyright © 2010-2024 three.js authors

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
```


## Film export (presentation films, docs/architecture/presentation-films.md)

`mp4-muxer.js` and `webm-muxer.js` are the packages' `build/*.mjs`
files, unmodified but renamed `.js` (every static server sends `.js`
as JavaScript; some send `.mjs` as a download, which a module import
refuses), from `https://registry.npmjs.org/mp4-muxer/-/mp4-muxer-5.2.2.tgz`
and `https://registry.npmjs.org/webm-muxer/-/webm-muxer-5.1.4.tgz`,
fetched 2026-10-15. `gifenc.esm.js` is `dist/gifenc.esm.js` from
`https://registry.npmjs.org/gifenc/-/gifenc-1.0.3.tgz` with its one
trailing `//# sourceMappingURL=gifenc.esm.js.map` comment removed (the
map is not shipped; upstream SHA-256
`217761244379253ba5815510d3048a4fb1f4c4dbe9fb51a5cfbfe26f34eec093`).
None imports anything or touches the network.

mp4-muxer:

```
MIT License

Copyright (c) 2023 Vanilagy

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```

webm-muxer:

```
MIT License

Copyright (c) 2022 Vanilagy

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```

gifenc:

```
The MIT License (MIT)
Copyright (c) 2017 Matt DesLauriers

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE
OR OTHER DEALINGS IN THE SOFTWARE.
```

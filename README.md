# popkern
Experimental BRep implementation built entirely by agents (in Rust)

## License

popkern is licensed under the **GNU Affero General Public License v3.0
only** (`LICENSE`) **with an additional permission** (`LICENSE-EXCEPTION`)
for user scripts and produced geometry.

### Licensing FAQ

**Do I have to publish scripts I write against popkern?** No. Part
scripts, Rhai/Python scripts, `.pkn` documents and migrations that use
popkern's published APIs are *User Scripts* under `LICENSE-EXCEPTION`:
they are not covered works, and they are excluded from the
Corresponding Source you would otherwise have to offer. License them
however you like.

**Who owns the models, meshes and STEP files popkern produces?** You
do. Output of any kind, produced locally or by a service you run, is
not a covered work and popkern's copyright holders claim no rights in
it.

**I run popkern behind a web service. What does AGPL §13 require?**
You must offer users of that service the Corresponding Source of the
*popkern* you run — including any changes you made to the kernel —
not your scripts, not your models.

**I modified the kernel / I ship a program that links popkern.** The
popkern portion stays AGPL: convey its source (with your kernel
modifications) as the license requires. Your own scripts are still
excluded by the exception. The exception is deliberately not a
Classpath-style linking exception — the kernel is copyleft; your work
on top of it is yours.

**Third-party dependencies.** Every crate in the dependency graph is
under a permissive license (see `THIRD-PARTY-NOTICES.md`, regenerated
by `cargo about generate`); `deny.toml` enforces the allowlist in CI.
See `docs/decisions/2026-08-30-user-scripts-and-output-license-exception.md`.

# Contributing to popkern

Thanks for wanting to help! popkern is early, and the most useful
things anyone can do right now are mostly not code.

## Break it

The single most valuable contribution is a part popkern can't build.
[Open an issue](https://github.com/popkern-dev/popkern/issues) with:

- **The script**, or the smallest version of it that still fails. A
  STEP file is great too, if it's an import that goes wrong.
- **What you expected**, and what happened instead. If popkern declined,
  paste the whole message, numbers and all. If it built something
  wrong without complaining, say so loudly: a quiet failure is the
  worst kind of bug, and I want to hear about every one.
- **Where you ran it**: popkern live in the browser, `popkern serve`,
  or the command line, plus the version (`popkern --version`) or the
  date you used the website.

Every failing part becomes a test case, so it stays fixed.

Found a security problem? Please don't open a public issue. See the
[security policy](https://popkern.dev/security/) for how to report it
privately.

## Code

Pull requests are welcome but not expected! For anything bigger than a
small fix, open an issue first so we can talk about the approach before
you spend time on it. Kernel changes in particular tend to have
consequences that aren't obvious from the outside.

[`AGENTS.md`](AGENTS.md) is the working guide for everyone changing
popkern, human or agent. It has the crate map, the dependency rules,
the code rules and the standard every change is held to. Read it
before you start.

Every change has to pass the same gates I run before every commit, and
CI checks them too:

```sh
cargo fmt --check
cargo clippy -- -D warnings
tools/quality-check.sh --full     # the kernel suite
tools/quality-check.sh --physics  # also, if you touched the analysis crates
```

The toolchain is pinned in `rust-toolchain.toml`, so rustup picks the
right one for you.

popkern is written by AI agents, so of course you're welcome to use
them too. The bar is the same either way: you're responsible for what
you submit, it has to pass the gates, and it has to make sense to a
reviewer.

## The contributor agreement

Before I can merge your first pull request, you'll need to sign the
[popkern Contributor License Agreement](CLA.md). The CLA Assistant bot
will ask you on your PR. It takes one comment, and you only do it
once.

Here's what it means in plain English, and why I'm asking.

**You keep your copyright.** The agreement is a licence, not a
transfer. It gives me the right to use your contribution in popkern
and to license it under any terms, including terms other than the
AGPL.

**Why I want that.** There are two futures I want to keep open. One is
relicensing popkern under something more permissive, like the MPL, if
that ever serves the people using it better. The other is selling
commercial licences to companies that can't use the AGPL, to pay for
popkern's development. Both are only possible if one person can
license all of the code. Once outside contributions are in without an
agreement like this one, that door closes for good.

**My promise in return.** Whatever else happens, popkern, including
your contribution, stays available as open source, with the exception
that leaves people's scripts and models theirs. A commercial licence
is only ever offered alongside the open source version, never instead
of it, and that promise is written into the agreement itself. Anything
already released under the AGPL stays that way forever, too.

If you're contributing as part of your job, check that your employer
is OK with it first. If you're contributing on behalf of a company,
get in touch before you open the PR: the company needs to sign its own
agreement.

## Licence

popkern is licensed under AGPL-3.0-only with the
[user scripts and output exception](LICENSE-EXCEPTION). Your
contribution is released under the same terms, plus the rights you
grant me in the [CLA](CLA.md).

## Say hi

Questions, ideas, or something you built? Come and say hi on the
[Discord](https://discord.gg/38p9xBBwff).

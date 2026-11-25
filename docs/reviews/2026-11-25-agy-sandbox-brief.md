# Brief: let Claude Code run `agy` (Antigravity CLI) for code reviews

## Goal
Claude Code sessions should be able to run `agy -p ...` non-interactively, for example for a
second-model code review alongside Claude's, using Brian's **Google account sign-in (Gemini
Pro subscription)**. Per-token API billing must not be used.

## What happened (2026-11-25, session in ~/Code/Julia/cooccurrence-analysis-r, auto mode)
1. `agy -p "Reply with the single word OK" --print-timeout 60s < /dev/null` ran inside the
   Bash sandbox. The network proxy refused outbound connections to:
   - `oauth2.googleapis.com:443` (sign-in / token refresh)
   - `antigravity-unleash.goog:443` (feature flags)
   - `play.googleapis.com:443` (telemetry; failures were only logged)
2. Claude re-ran the command with
   `allowed_domains: ["oauth2.googleapis.com", "antigravity-unleash.goog", "play.googleapis.com", "*.googleapis.com"]`.
   The **auto-mode permission check denied it as "Unauthorized Access."** The likely reason is
   that `agy` would use stored Google credentials that the session had not been given.
3. The command never got far enough to show whether `agy` was signed in.

So there are **two separate barriers**: (a) the sandbox's network allowlist, and (b) the
auto-mode permission check's judgment about using account credentials. A fix for (a) alone may
not get past (b).

## Constraints
- **Never use the `gemini` CLI for this.** `~/.gemini/settings.json` has
  `security.auth.selectedType: "gemini-api-key"`, which bills per token. Do not set
  `GEMINI_API_KEY`, and do not change that setting as a "fix."
- `agy` signs in with the Google account (OAuth). Keep it that way.
- Keep the change narrow: only the hosts `agy` needs, and only the `agy` command. Prefer an
  `agy`-specific permission over loosening the sandbox in general.
- Reviews should be read-only. Claude invokes `agy` with `--mode plan`.

## What to implement (in Claude Code settings)
1. **Network:** add the hosts `agy` needs to the sandbox's network allowlist. Start with the
   three above. The model-serving host(s) are **not yet known**: the test never got past
   sign-in. Find them by running the acceptance test below and reading the
   `<sandbox_violations>` block. Do not guess them.
2. **Permission:** add an allow rule for the `agy` command (for example `Bash(agy:*)`, or a
   narrower `Bash(agy -p:*)`), so auto mode treats it as authorized rather than as
   unauthorized credential use.
3. **Alternative**, if the two above aren't enough: exclude `agy` from the sandbox, the way
   `gh` already is (managed settings say `gh` runs as a standalone command via
   `excludedCommands` so it can reach the Keychain). This is broader; use it only if needed.
4. Choose the scope deliberately: the user-level `~/.claude/settings.json` (all projects) or
   per-project settings.

## Acceptance test
From a fresh Claude Code session in auto mode, Claude runs:

    agy -p "Reply with the single word OK" --mode plan --print-timeout 90s < /dev/null

Passes when: the output contains `OK`, there is no `<sandbox_violations>` block other than
optional `play.googleapis.com` telemetry, the auto-mode check does not deny it, and `gemini`
or `GEMINI_API_KEY` is not involved. If `agy` reports it is signed out, Brian runs `agy` once
in an ordinary terminal to sign in, and then the test is repeated.

## Related (optional)
The `mmreview` skill (`~/.claude/skills/mmreview/SKILL.md`) calls the `gemini` CLI, which is
in API-key mode, and wraps commands in `timeout`, which macOS lacks (no coreutils
`gtimeout` installed). Its 120 s limit is also too short for large diffs. If `mmreview` is to
be the standard multi-model review, change it to call `agy -p ... --mode plan --output-format json`
in place of `gemini`, and replace `timeout` with `agy`'s own `--print-timeout`.

#[derive(Debug, Clone, PartialEq, Eq)]
pub(super) enum RunnerCommand {
    Run(RunnerSelection),
    ListFixtures,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub(super) struct RunnerSelection {
    pub preset: FixturePreset,
    pub fixture_stems: Vec<String>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum FixturePreset {
    Smoke,
    Full,
}

impl FixturePreset {
    pub(super) fn label(self) -> &'static str {
        match self {
            Self::Smoke => "smoke",
            Self::Full => "full",
        }
    }
}

pub(super) const SMOKE_FIXTURES: &[&str] = &["t0", "tee", "cube"];

pub(super) fn parse_runner_command<I>(args: I) -> Result<RunnerCommand, String>
where
    I: IntoIterator<Item = String>,
{
    let mut preset = FixturePreset::Smoke;
    let mut fixture_stems = Vec::new();
    let mut iter = args.into_iter().skip(1);
    while let Some(arg) = iter.next() {
        match arg.as_str() {
            "--preset" => {
                let value = iter
                    .next()
                    .ok_or_else(|| "--preset requires smoke or full".to_string())?;
                preset = parse_preset(&value)?;
            }
            "--fixture" => {
                let value = iter
                    .next()
                    .ok_or_else(|| "--fixture requires a fixture stem".to_string())?;
                push_fixture_stems(&mut fixture_stems, &value);
            }
            "--fixtures" => {
                let value = iter.next().ok_or_else(|| {
                    "--fixtures requires comma-separated fixture stems".to_string()
                })?;
                push_fixture_stems(&mut fixture_stems, &value);
            }
            "--list-fixtures" => return Ok(RunnerCommand::ListFixtures),
            "--help" | "-h" => return Err(usage()),
            other => {
                return Err(format!("unknown argument `{other}`\n{}", usage()));
            }
        }
    }
    Ok(RunnerCommand::Run(RunnerSelection {
        preset,
        fixture_stems,
    }))
}

pub(super) fn selected_fixture_stems<'a>(
    selection: &'a RunnerSelection,
    all_stems: &'a [&'static str],
) -> Result<Vec<&'a str>, String> {
    if selection.fixture_stems.is_empty() {
        return match selection.preset {
            FixturePreset::Smoke => Ok(SMOKE_FIXTURES.to_vec()),
            FixturePreset::Full => Ok(all_stems.to_vec()),
        };
    }
    let mut selected = Vec::with_capacity(selection.fixture_stems.len());
    for stem in &selection.fixture_stems {
        if !all_stems.contains(&stem.as_str()) {
            return Err(format!(
                "unknown fixture `{stem}`\n{}",
                fixture_list(all_stems)
            ));
        }
        if !selected.contains(&stem.as_str()) {
            selected.push(stem.as_str());
        }
    }
    Ok(selected)
}

pub(super) fn fixture_list(all_stems: &[&'static str]) -> String {
    format!("available fixtures: {}", all_stems.join(", "))
}

fn parse_preset(value: &str) -> Result<FixturePreset, String> {
    match value {
        "smoke" => Ok(FixturePreset::Smoke),
        "full" => Ok(FixturePreset::Full),
        _ => Err(format!("unknown preset `{value}`; expected smoke or full")),
    }
}

fn push_fixture_stems(out: &mut Vec<String>, value: &str) {
    out.extend(
        value
            .split(',')
            .map(str::trim)
            .filter(|stem| !stem.is_empty())
            .map(str::to_string),
    );
}

fn usage() -> String {
    [
        "usage: pns_toms_step_examples [--preset smoke|full]",
        "       pns_toms_step_examples --fixture tee [--fixture cube]",
        "       pns_toms_step_examples --fixtures tee,cube",
        "       pns_toms_step_examples --list-fixtures",
    ]
    .join("\n")
}

#[cfg(test)]
mod tests {
    use super::*;

    const ALL: &[&str] = &["t0", "tee", "cube", "suzanne_all_configurations"];

    #[test]
    fn default_selection_is_bounded_smoke_set() {
        let RunnerCommand::Run(selection) =
            parse_runner_command(["bin".to_string()]).expect("default args")
        else {
            panic!("default command must run");
        };
        assert_eq!(
            selected_fixture_stems(&selection, ALL).unwrap(),
            vec!["t0", "tee", "cube"]
        );
    }

    #[test]
    fn explicit_fixtures_override_preset_and_deduplicate() {
        let RunnerCommand::Run(selection) = parse_runner_command([
            "bin".to_string(),
            "--preset".to_string(),
            "full".to_string(),
            "--fixtures".to_string(),
            "tee,cube,tee".to_string(),
        ])
        .expect("explicit fixtures") else {
            panic!("fixture command must run");
        };
        assert_eq!(
            selected_fixture_stems(&selection, ALL).unwrap(),
            vec!["tee", "cube"]
        );
    }
}

# Linting

Two linters check the code in the GitLab CI pipeline:
[Fortitude](https://github.com/PlasmaFAIR/fortitude) for the Fortran sources in `src` and
[Ruff](https://docs.astral.sh/ruff/) for the Python scripts in `docs` and `tools`.
A pipeline fails when the `fortitude:` or `ruff:` stage reports a violation.

## Fortitude (Fortran)

[Fortitude is a Fortran linter](https://github.com/PlasmaFAIR/fortitude) that can be installed from PYPI via

    # With uv:
    uv tool install fortitude-lint@latest

    # With pip:
    pip install fortitude-lint

After installing, check the version

    fortitude --version

Fortitude is used by navigating to the piclas repository and running

    cd ~/piclas
    fortitude check --output-format=grouped src

which lists any errors that are encountered with their corresponding violation ID.

The Fortitude GitLab CI/CD stage `fortitude:` is defined in `.gitlab-ci.yml` and the Fortitude settings are stored in `fpm.toml`.

### Fortitude pre-commit hook
TODO: The pre-commit hook for Fortitude is configured in ...

### Fixing Fortitude violations

When a pipeline fails in the `fortitude` stage, an error is displayed showing the ID of the violation

    src/timedisc/timedisc_TimeStepECIM.f90:
    606:42 C141 'exit' statement in named 'do' loop missing label 'SUBROUTINE ExactPushSingleParticle(iPart, dt)

which in this case is "C141" and details on this violation can be found by running

    fortitude explain C141

which returns

    C141: missing-exit-or-cycle-label

    Fix is sometimes available.

    What does it do?
    When using exit or cycle in a named do loop, the exit/cycle statement
    should use the loop name

    Example
    name: do
      exit name
    end do name

    Using named loops is particularly useful for nested or complicated loops, as it
    helps the reader keep track of the flow of logic. It's also the only way to exit
    or cycle outer loops from within inner ones.

Fortitude can automatically fix many linter warnings and errors. Simply navigate to the piclas repository and use the `--fix` flag
for Fortitude

    cd ~/piclas
    fortitude check --output-format=grouped src --fix

or

    cd ~/piclas
    fortitude check --output-format=grouped src/particles/pic/models/pic_models.f90 --select C141 --fix

to only fix specific violation IDs in specific .f90 files.

Note that there are two flags for auto-fixing `--fix` and `--unsafe-fixes`
| Flag                       | Description                                                                                |
| :------------------------- | :----------------------------------------------------------------------------------------- |
| `check`                    | Only checks for issues without applying any fixes.                                         |
| `--fix`                    | Applies **safe** automatic fixes (e.g., style-based rules).                                |
| `--unsafe-fixes`           | Applies **both safe and unsafe** automatic fixes (use with caution).                       |
| `--diff`                   | Shows the changes that would be made by `--fix` or `--unsafe-fixes` without applying them. |
| `--help`                   | Displays help information about Fortitude commands and flags.                              |
| `--version`                | Shows the installed version of Fortitude.                                                  |
| `--output-format=<format>` | Specifies the output format (e.g., `grouped`, `json`, `sarif`).                            |

Unsafe fixes may change program behavior or introduce new issues, so use with caution.

## Ruff (Python)

[Ruff is a Python linter and code formatter](https://docs.astral.sh/ruff/) that can be installed from PYPI via

    # With uv:
    uv tool install ruff

    # With pip:
    pip install ruff

After installing, check the version

    ruff --version

Ruff is used by navigating to the piclas repository and running

    cd ~/piclas
    ADDIGNORE=,F821,PLW0602,UP031,DTZ011,TRY002,BLE001,E722,SIM115,RUF012,SIM102,B018,EXE001,S110,LOG015,TRY300,PLC0206,PLR1704,F509,PLE1300,PLW0127,S102,SIM113,RET504,F841,C408,FLY002
    ruff check --extend-ignore=E201,E202,E203,E221,E222,E225,E231,E271,E272"$ADDIGNORE" --line-length=132 --preview docs tools --statistics

which lists any errors that are encountered with their corresponding violation ID, for example

    E722 Do not use bare `except`
      --> tools/visuBox_PartState/visuBox_PartState.py:59:1
       |
    57 |     time = float(timestr)
    58 |     newFile = re.sub(timestr+'.h5', '', args.statefile)+timestr+'_visuBox.h5'
    59 | except :
       | ^^^^^^
    60 |     print("something wrong with "+args.statefile)
    61 |     sys.exit(0)
       |

Running two command with `--statistics` will list a summary with the number of hits for each violation and with a short description

    ruff check --extend-ignore=E201,E202,E203,E221,E222,E225,E231,E271,E272"$ADDIGNORE" --line-length=132 --preview docs tools --statistics

which gives, for example, the following

    207     F821    undefined-name
     62     UP031   printf-string-formatting
     50     PLW0602 global-variable-not-assigned
     24     TRY002  raise-vanilla-class
     21     DTZ011  call-date-today
     19     E722    bare-except
     19     SIM115  open-file-with-context-handler
     17     BLE001  blind-except
      8     SIM102  collapsible-if
      7     RUF012  mutable-class-default
      6     F841    unused-variable
      5     B018    useless-expression
      5     EXE001  shebang-not-executable
      5     S110    try-except-pass

The Ruff GitLab CI/CD stage `ruff:` is defined in `.gitlab-ci.yml`.

### Ruff pre-commit hook
TODO: The pre-commit hook for Ruff is configured in ...

### Fixing Ruff violations

When a pipeline fails in the `ruff` stage, an error is displayed showing the ID of the violation

    FLY002 Consider f-string instead of string join
      --> tools/paraview/pdf/43/python_filter_generator.py:81:12
       |
    79 |       </StringVectorProperty>''' % requestUpdateExtent
    80 |
    81 |     return '\n'.join([requestData, requestInformation, requestUpdateExtent])
       |            ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
       |
    help: Replace with f-string

    Found 1 error.

which in this case is "FLY002" and details on this violation can be found by running

    ruff rule FLY002

which returns

    # static-join-to-f-string (FLY002)

    Derived from the **flynt** linter.

    Fix is always available.

    ## What it does
    Checks for `str.join` calls that can be replaced with f-strings.

    ## Why is this bad?
    f-strings are more readable and generally preferred over `str.join` calls.

    ## Example
    ```python
    " ".join((foo, bar))
    ```

    Use instead:
    ```python
    f"{foo} {bar}"
    ```

    ## Fix safety
    The fix is always marked unsafe because the evaluation of the f-string
    expressions will default to calling the `__format__` method of each
    object, whereas `str.join` expects each object to be an instance of
    `str` and uses the corresponding string. Therefore it is possible for
    the values of the resulting strings to differ, or for one expression
    to raise an exception while the other does not.

    ## References
    - [Python documentation: f-strings](https://docs.python.org/3/reference/lexical_analysis.html#f-strings)

Ruff can automatically fix many violations by adding the `--fix` flag

    cd ~/piclas
    ruff check --extend-ignore=E201,E202,E203,E221,E222,E225,E231,E271,E272"$ADDIGNORE" --line-length=132 --preview docs tools --fix

which might return with

    No fixes available (1 hidden fix can be enabled with the `--unsafe-fixes` option).

and in this case an unsafe fix might work via `--unsafe-fixes`

    cd ~/piclas
    ruff check --extend-ignore=E201,E202,E203,E221,E222,E225,E231,E271,E272"$ADDIGNORE" --line-length=132 --preview docs tools --unsafe-fixes

which shows a diff of how the violation would be fixed

```diff
78 |         extent that your filter ask up stream for.</Documentation>
79 |       </StringVectorProperty>''' % requestUpdateExtent
80 |
   -     return '\n'.join([requestData, requestInformation, requestUpdateExtent])
81 +     return f'{requestData}\n{requestInformation}\n{requestUpdateExtent}'
82 |
83 |
84 |
note: This is an unsafe fix and may change runtime behavior

Found 1 error.
[*] 1 fixable with the `--fix` option
```

and can be applied by running the command with both flags `--fix --unsafe-fixes`

    ruff check --extend-ignore=E201,E202,E203,E221,E222,E225,E231,E271,E272"$ADDIGNORE" --line-length=132 --preview docs tools --fix --unsafe-fixes

It is important to make sure that the code still behaves as intended when using an unsafe fix.
| Flag | Description |
|------|-------------|
| `--fix` | Automatically fix lint errors where possible. |
| `--watch` | Watch for file changes and re-run Ruff automatically. |
| `--select=<CODE>` | Select specific lint rules to apply (e.g., `E`, `F`, `I`, `UP`). |
| `--ignore=<CODE>` | Ignore specific lint rules. |
| `--line-length=<N>` | Set the maximum line length (default: `88`). |
| `--statistics` | Show statistics about the linting process. |
| `--diff` | Show a diff of changes Ruff would make (useful for CI/CD). |

## Finding earlier fixes

If the fix is not straightforward, look at how other developers fixed the same violation ID.
Search the complete git history (commit messages and changes to the code) via

    cd ~/piclas
    git log --all --grep="C141" -p

and search the source code for places where the violation is suppressed via

    grep -rin --include=*.f90 C141 src
    grep -rin --include=*.py S112 docs tools

A violation is either fixed in the code or, for false positives, suppressed with a comment:

- **Fortitude**: add `! allow(C141)` in the line before the one with the trigger.
  False positives are common, because Fortitude still has shortcomings regarding macros and pre-processor statements.

  ```diff
     IF(PDM%ParticleInside(iPart)) THEN
       ASSOCIATE ( oldSpec => PartSpecies(iPart) ,&
             newSpec => SpecDSMC(PartSpecies(iPart))%NextIonizationSpecies )
  +      ! allow(C141) because of false-positive Fortitude check
         IF(newSpec.EQ.0) CYCLE
  ```

- **Ruff**: add `# noqa: S112` at the end of the line with the trigger.

      except :  # noqa: S112
          continue

## Separate commits

Fix linter violations in a separate commit and name the violation IDs (e.g. C141 or FLY002) in the commit message,
so that other developers find the fix with `git log --grep` when they tackle the same issue.

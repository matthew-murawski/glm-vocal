# Agent Instructions

Consult BLUEPRINT.md for development blueprint, or SPEC.md!

## Commenting in MATLAB Files
- keep comments in lowercase, with uppercase allowed only for abbreviations (e.g., `ACC`).
- provide two types of comments:
  - section-by-section comments that explain the intention of the section, and walk the user through what is being done. these will often be full sentences across multiple lines. do not prepend markers like "section:".
  - in-section comments (keep more concise)
- avoid obvious or trivial comments.
- do not over-comment.
- keep comments informal yet helpful.

Use this exact command to run the test suite:
./scripts/run_matlab_tests.sh
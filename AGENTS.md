# Agent Instructions

- Users will always have the structure "P" (containing paths) in their enviornment. Never override this structure.

- First, determine if you are writing a script or a function. A script is a one-click file that sets parameters (which are modifyable by the user) and calls an orchestrator function to do its bidding. Usually the parameters are monkey name and session number or something like that. A function is, well, a function that takes parameters and always has "function" in the first line. If it is a script, go to " Scripts" section. If it is a function, go to the "Functions" section of this doc.

## Scripts

Scripts must ALWAYS take the stereotyped form found in docs/script_template.txt, pasted here:

%% [script_file_name]
% [one-line description of the script]

clr; clc;

%% User Parameters
[Parameter = parameter;]
[OtherParameter = 3;]

[ListOfParameters = {'Param1', 'Param2', 'Param3'};]

%% Execute Analysis
[fprintf(['Starting whatever for %s | session %d]);]

[the actual runner script being called]

[fprintf('Whatever has been generated or whatever. Done.\n');]

items found in [] can be modified. Items not encased in [] must always be the same for each script.

Keep comments on scripts light.

You do not need to test scripts.

Scripts must ALWAYS start with "run_". They must always be placed in the scripts folder. If something is placed in the scripts folder, it must function as a one-click runner.

## Functions

- keep comments in lowercase, with uppercase allowed only for abbreviations (e.g., `ACC`).
- structure every function into clear sections, each preceded by a section comment (multiple sentences if necessary) that walks the reader through the next block of logic. most functions will have several sections; very small helpers may only need one.
- it is quite alright if a section comment requires multiple lines. really make sure you tell the user about why you're writing a section if it's not totally clear. this usually requires multiple lines of a comment (in sentence form) to get right. this will not be more than 1/3 of all comments though.
- add concise in-section comments when a line or two needs extra context, but avoid restating the code.
- do not over-comment and skip obvious remarks; focus on intent, preconditions, and tricky details.
- keep comments informal yet helpful.

## Testing
Use this exact command to run the test suite:
./scripts/run_matlab_tests.sh

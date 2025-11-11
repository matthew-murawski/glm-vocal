# Agent Instructions

Hi. You are helping the user write code in MATLAB for a marmoset vocalization analysis project.

First, check SPEC.md to see what you are building. 
Next, check CHECKLIST.md to see where you are at in the process.
As you go, check off items in CHECKLIST.md by changing [ ] to [X]

Use this exact command to run the test suite:
./scripts/run_matlab_tests.sh

When you finish implementing a prompt, make sure to commit to git.

# Information
- Users will always have the structure "P" (containing paths) in their enviornment. Never override this structure.

- First, determine if you are writing production code (or what will BECOME production code) or dev code. Dev code includes human-verification scripts (that you write for each "phase" for the human to verify you have correctly implemented the prompt) as well as MATLAB unit/e2e tests. Production code is either a script or a function that will be used in the final product.

If it is dev code, please still write good, principled code; however, you don't need to strictly follow the below rules (especially you don't need to follow commenting rules). However, if it is production code, please follow the following rules.

# Production Code

If it is production code, first, determine if you are writing a script or a function. A script is a one-click file that sets parameters (which are modifyable by the user) and calls an orchestrator function to do its bidding. Usually the parameters are monkey name and session number or something like that. A function is, well, a function that takes parameters and always has "function" in the first line. If it is a production script, go to " Scripts" section. If it is a function, go to the "Functions" section of this doc.

## Scripts

Production scripts must ALWAYS take the stereotyped form found in docs/script_template.txt, pasted here:

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

Whenever you add a new script, you need to register it to the registry. To do so, you need to run this function (modifying the parameters accordingly):

success = register_script(...
    'name', 'Single Session Firing Rate', ...
    'id', 'firing_rate_session', ...
    'category', 'Spike', ...
    'subcategory', 'FiringRate', ...
    'script_path', 'scripts/runners/spike/firing_rate_session.m', ...
    'description', 'Analyze firing rates for a single session', ...
    'order', 10);

the output, "success", will be true if registration is successful, and false otherwise.

## Functions

- keep comments in lowercase, with uppercase allowed only for abbreviations (e.g., `ACC`).
- structure every function into clear sections, each preceded by a section comment (multiple sentences if necessary) that walks the reader through the next block of logic. most functions will have several sections; very small helpers may only need one.
- it is quite alright if a section comment requires multiple lines. really make sure you tell the user about why you're writing a section if it's not totally clear. this usually requires multiple lines of a comment (in sentence form) to get right. this will not be more than 1/3 of all comments though.
- add concise in-section comments when a line or two needs extra context, but avoid restating the code.
- do not over-comment and skip obvious remarks; focus on intent, preconditions, and tricky details.
- keep comments informal yet helpful.
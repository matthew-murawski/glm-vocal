function tests = test_struct_tools
%TEST_STRUCT_TOOLS Unit tests for struct utility helpers.
%   tests = TEST_STRUCT_TOOLS() returns function-based tests.

% register local test functions with the MATLAB framework
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
% set up reusable handles so we only instantiate helpers once
testCase.TestData.tools = struct_tools();
end

function testStructMergeDeep(testCase)
% grab helper bundle and prepare nested structs to merge
tools = testCase.TestData.tools;
base = struct('a', 1, 'b', struct('c', 2, 'd', 3));
override = struct('b', struct('d', 5, 'e', 6), 'f', 7);

% perform the merge and capture the composite struct
merged = tools.merge(base, override);
expectedB = struct('c', 2, 'd', 5, 'e', 6);

% ensure overrides land where expected while preserving other fields
assert(isequal(merged.a, 1));
assert(isequal(merged.b, expectedB));
assert(isequal(merged.f, 7));
end

function testStructMergeNonStructOverride(testCase)
% set up a base struct and an overriding scalar value
tools = testCase.TestData.tools;
base = struct('a', struct('nested', 1));
override = struct('a', 10);

% confirm non-struct overrides replace the branch entirely
merged = tools.merge(base, override);
assert(isequal(merged.a, 10));
end

function testStructGetWithDefault(testCase)
% design a nested struct and query two different paths
tools = testCase.TestData.tools;
data = struct('alpha', struct('beta', struct('gamma', 42)));

% fetch a valid path and a missing path to exercise defaults
value = tools.get(data, 'alpha.beta.gamma', NaN);
missing = tools.get(data, 'alpha.beta.delta', 'missing');
assert(isequal(value, 42));
assert(strcmp(missing, 'missing'));
end

function testStructSetNested(testCase)
% start from an empty struct and set a deep path
tools = testCase.TestData.tools;
data = struct();

% ensure intermediate structs appear and value sticks
data = tools.set(data, 'outer.inner.value', 11);
assert(isfield(data, 'outer'));
assert(isfield(data.outer, 'inner'));
assert(isfield(data.outer.inner, 'value'));
assert(isequal(data.outer.inner.value, 11));
end

function testStructSetOverwriteNonStruct(testCase)
% begin with a scalar field to verify overwrite behaviour
tools = testCase.TestData.tools;
data = struct('outer', 5);

% applying set should convert the field into a struct container
data = tools.set(data, 'outer.inner', 12);
assert(isstruct(data.outer));
assert(isequal(data.outer.inner, 12));
end

function testSplitPathValidation(testCase)
% pass an invalid path type and expect a descriptive failure
tools = testCase.TestData.tools;
verifyError(testCase, @() tools.get(struct(), 123, []), 'struct_tools:InvalidPath');
end

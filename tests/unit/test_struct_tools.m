function tests = test_struct_tools
%TEST_STRUCT_TOOLS Unit tests for struct utility helpers.
%   tests = TEST_STRUCT_TOOLS() returns function-based tests.

tests = functiontests(localfunctions);
end

function setupOnce(testCase)
%SETUPONCE Cache helper handles for convenience.
testCase.TestData.tools = struct_tools();
end

function testStructMergeDeep(testCase)
tools = testCase.TestData.tools;
base = struct('a', 1, 'b', struct('c', 2, 'd', 3));
override = struct('b', struct('d', 5, 'e', 6), 'f', 7);
merged = tools.merge(base, override);
expectedB = struct('c', 2, 'd', 5, 'e', 6);
assert(isequal(merged.a, 1));
assert(isequal(merged.b, expectedB));
assert(isequal(merged.f, 7));
end

function testStructMergeNonStructOverride(testCase)
tools = testCase.TestData.tools;
base = struct('a', struct('nested', 1));
override = struct('a', 10);
merged = tools.merge(base, override);
assert(isequal(merged.a, 10));
end

function testStructGetWithDefault(testCase)
tools = testCase.TestData.tools;
data = struct('alpha', struct('beta', struct('gamma', 42)));
value = tools.get(data, 'alpha.beta.gamma', NaN);
missing = tools.get(data, 'alpha.beta.delta', 'missing');
assert(isequal(value, 42));
assert(strcmp(missing, 'missing'));
end

function testStructSetNested(testCase)
tools = testCase.TestData.tools;
data = struct();
data = tools.set(data, 'outer.inner.value', 11);
assert(isfield(data, 'outer'));
assert(isfield(data.outer, 'inner'));
assert(isfield(data.outer.inner, 'value'));
assert(isequal(data.outer.inner.value, 11));
end

function testStructSetOverwriteNonStruct(testCase)
tools = testCase.TestData.tools;
data = struct('outer', 5);
data = tools.set(data, 'outer.inner', 12);
assert(isstruct(data.outer));
assert(isequal(data.outer.inner, 12));
end

function testSplitPathValidation(testCase)
tools = testCase.TestData.tools;
verifyError(testCase, @() tools.get(struct(), 123, []), 'struct_tools:InvalidPath');
end

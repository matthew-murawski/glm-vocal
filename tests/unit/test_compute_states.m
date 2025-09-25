function tests = test_compute_states
% exercise compute_states against response-required conversational logic.

%% register the test cases
% we expose the local functions to matlab's function-based test harness.
tests = functiontests(localfunctions);
end

function stim = makeStim(dt, nBins)
% build a reusable stimulus timeline for the tests.

%% define the timeline grid
% we create an evenly spaced grid with a fully valid mask by default.
t = (0:nBins-1)' * dt;
stim = struct('t', t, 'dt', dt, 'mask', struct('good', true(nBins, 1)));
end

function cfg = makeStateCfg()
% helper to generate a standard state configuration.

%% expose response and state windows aligned with the session rules
cfg = struct('response_window_s', 5.0, 'state_window_s', 5.0);
end

function ev = makeEvent(kind, tOn, tOff)
% convenience constructor for event structs used in tests.

%% pack the supplied values with an empty label to mirror repo expectations
ev = struct('kind', kind, 't_on', tOn, 't_off', tOff, 'label', "");
end

function testRequiresProducedResponse(testCase)
% heard calls without a timely produced response should remain spontaneous.

%% craft a perceived event with a late produced call beyond the response window
stim = makeStim(1.0, 15);
ev = [ ...
    makeEvent('perceived', 1.0, 1.3); ...
    makeEvent('produced', 7.0, 7.3) ...
];

%% compute states and ensure no conversational bins are marked
states = compute_states(ev, stim, makeStateCfg());
expectedConvo = false(15, 1);
expectedSpon = true(15, 1);
verifyEqual(testCase, states.convo, expectedConvo);
verifyEqual(testCase, states.spon, expectedSpon);
end

function testConversationChainsWithReplies(testCase)
% alternating reply sequences within the window sustain the conversation.

%% assemble interleaved perceived/produced sequences obeying the 5s rule
stim = makeStim(1.0, 16);
ev = [ ...
    makeEvent('perceived', 2.0, 2.2); ...
    makeEvent('produced', 4.0, 4.3); ...
    makeEvent('produced', 5.0, 5.2); ...
    makeEvent('perceived', 7.0, 7.4); ...
    makeEvent('produced', 9.0, 9.2) ...
];

%% compute states and confirm the conversation spans the expected bins
states = compute_states(ev, stim, makeStateCfg());
expectedConvo = false(16, 1);
expectedConvo(3:15) = true; % conversation from 2s through 15s (5s post reply)
verifyEqual(testCase, states.convo, expectedConvo);
verifyEqual(testCase, states.spon, ~expectedConvo);
end

function testConversationEndsWithoutReplyToHeard(testCase)
% missing a produced response to a heard reply should end the conversation.

%% build a conversation that fails to respond to the next heard sequence
stim = makeStim(1.0, 18);
ev = [ ...
    makeEvent('perceived', 1.0, 1.2); ...
    makeEvent('produced', 3.0, 3.2); ...
    makeEvent('perceived', 6.0, 6.3) ...
];

%% compute states and verify the window closes five seconds after the last heard start
states = compute_states(ev, stim, makeStateCfg());
expectedConvo = false(18, 1);
expectedConvo(2:11) = true; % from 1s through 11s (5s after last heard start)
verifyEqual(testCase, states.convo, expectedConvo);
verifyEqual(testCase, states.spon, ~expectedConvo);
end

function testConversationEndsWithoutNextHeard(testCase)
% lack of a follow-up heard call terminates the conversation after 5s.

%% create a single heard→produced exchange with no subsequent heard reply
stim = makeStim(1.0, 14);
ev = [ ...
    makeEvent('perceived', 0.5, 0.7); ...
    makeEvent('produced', 2.0, 2.3) ...
];

%% compute states and ensure the conversation extends 5s past the last produced call
states = compute_states(ev, stim, makeStateCfg());
expectedConvo = false(14, 1);
expectedConvo(1:8) = true; % from 0.5s through 7.5s (5s post produced)
verifyEqual(testCase, states.convo, expectedConvo);
verifyEqual(testCase, states.spon, ~expectedConvo);
end

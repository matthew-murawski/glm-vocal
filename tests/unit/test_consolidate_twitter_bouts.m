function tests = test_consolidate_twitter_bouts
% section registration
% expose consolidation helper tests to the matlab runner.
tests = functiontests(localfunctions);
end

function testNoTwitterReturnsInput(testCase)
% section passthrough
% verify non-twitter events are not altered.
events = makeEvents([0, 0.5], {'produced', 'perceived'}, {'phee', 'heard'});
out = consolidate_twitter_bouts(events, 1.5);

testCase.verifyEqual(out, events);
end

function testTwitterWithinWindowMerges(testCase)
% section merge window
% ensure multiple twitter syllables within the window collapse to one event.
events = makeEvents([0, 0.4, 2.0], {'produced', 'produced', 'produced'}, {'twitterA', 'twitterB', 'twitterC'});
out = consolidate_twitter_bouts(events, 1.0);

testCase.verifyEqual(numel(out), 2);
ons = [out.t_on];
testCase.verifyEqual(ons, [0, 2.0], 'AbsTol', 1e-12);
end

function testTwitterBeyondWindowStartsNewBout(testCase)
% section new bout
% confirm onsets separated by more than the window create a new event.
events = makeEvents([0, 1.6, 3.3], {'produced', 'produced', 'produced'}, {'twitter', 'twitter', 'twitter'});
out = consolidate_twitter_bouts(events, 1.5);

testCase.verifyEqual(numel(out), 3);
ons = [out.t_on];
testCase.verifyEqual(ons, [0, 1.6, 3.3], 'AbsTol', 1e-12);
end

function events = makeEvents(times, kinds, labels)
% helper to build a struct array of events
n = numel(times);
events = repmat(struct('kind', '', 't_on', 0, 't_off', 0, 'label', ""), n, 1);
for ii = 1:n
    events(ii).kind = kinds{ii};
    events(ii).t_on = times(ii);
    events(ii).t_off = times(ii) + 0.1;
    events(ii).label = labels{ii};
end
end

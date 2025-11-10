function rec = empty_event_record()
% provide a consistent empty record template for event structures
rec = struct('kind', '', 't_on', 0, 't_off', 0, 'label', "");
end
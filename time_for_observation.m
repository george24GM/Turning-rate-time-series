function result= time_for_observation(n, time2)
    % n is the index of the observation in EEG for which you want to know the
    % corresponding time. Example: EEG is array of recordings, and you want to know
    % the time corresponding to the 13400000th observation. Then n=13400000.
    % time2 = 'hh:mm:ss' is the starting time of the recordings (find it in the
    % recording.txt reader).
    
    % Given parameters
    sampling_rate = 512; % Observations per second
    
    % Calculate total seconds elapsed
    total_seconds = n / sampling_rate;

    % Convert total seconds to hours, minutes, and seconds
    h1 = floor(total_seconds / 3600);
    carry_seconds = mod(total_seconds, 3600);
    m1 = floor(carry_seconds / 60);
    s1 = mod(carry_seconds, 60);

    % Splitting starting time into hours, minutes, and seconds
    [h2, m2, s2] = splitTime(time2);

    % Add seconds
    s = mod(s1 + s2, 60);
    carry_minutes = floor((s1 + s2) / 60);

    % Add minutes
    m = mod(m1 + m2 + carry_minutes, 60);
    carry_hours = floor((m1 + m2 + carry_minutes) / 60);

    % Add hours
    h = mod(h1 + h2 + carry_hours, 24);

    % Formatting the result
    result = sprintf('%02d:%02d:%02d', h, m, floor(s));
    disp(result);
end

function [hours, minutes, seconds] = splitTime(time)
    % Splitting time string into hours, minutes, and seconds
    parts = strsplit(time, ':');
    hours = str2double(parts{1});
    minutes = str2double(parts{2});
    seconds = str2double(parts{3});
end

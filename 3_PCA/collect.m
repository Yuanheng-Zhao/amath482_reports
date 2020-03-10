% Collect and arrange the data from three cameras into a single matrix
% Thus the matrix incorporates data points of all measurement types 
% taken by each camera over time.
function collected_data = collect(data1, data2, data3)
    [M,I] = min(data1(1:20,2));
    data1 = data1(I:end,:);
    [M,I] = min(data2(1:20,2));
    data2 = data2(I:end,:);
    [M,I] = min(data3(1:20,2));
    data3 = data3(I:end,:);

    % Trim the data to make them a consistent length.
    % For Test1, 2, and 3, the video recorded by camera 1 is (always) the shortest. 
    % Thus trim the other two as the length of video 1.
%     disp(length(data1));
%     disp(length(data2));
%     disp(length(data3));
    data2 = data2(1:length(data1), :);
    data3 = data3(1:length(data1), :);

    collected_data = [data1'; data2'; data3'];
end

function output_txt = datacursor(~,event_obj,labels)
% Display the position of the data cursor
% obj          Currently not used (empty)
% event_obj    Handle to event object
% output_txt   Data cursor text string (string or cell array of strings).
    pos = get(event_obj,'Position');
    idx = get(event_obj, 'DataIndex');
%     output_txt = {['X: ',num2str(pos(1),4)],...
%                  ['Y: ',num2str(pos(2),4)] };
%     output_txt{end+1} = ['Index: ', num2str(idx)];
    output_txt = [['Val: ',num2str(pos(2),4)],['Voxels: ', labels(idx)]];
end
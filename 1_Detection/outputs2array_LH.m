function arrout = outputs2array(fn, ixsOutputs)
    output_cell = cell(1,max(ixsOutputs)); 
    [output_cell{:}] = (fn());             
    arrout = output_cell(ixsOutputs);
    arrout = [arrout{:}];
end

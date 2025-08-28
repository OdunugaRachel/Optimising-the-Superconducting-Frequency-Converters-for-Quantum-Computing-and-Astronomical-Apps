function fcalc = removeHarmonics(fcalc, pumps, Idc)

    pumps3WM = pumps(:)+pumps(:).';
    pumps4WM = pumps3WM(:)+pumps(:).';

    if Idc == 0
        exclude = cat(2, pumps, transpose(unique(pumps4WM)));
    else
        exclude = unique(cat(2, transpose(unique(pumps3WM)), transpose(unique(pumps4WM)), pumps, pumps/2));
    end

    for i=1:length(fcalc)
        if ismember(fcalc(i),exclude)
            fcalc(i) = fcalc(i) + 1e3;
        end
    end

end
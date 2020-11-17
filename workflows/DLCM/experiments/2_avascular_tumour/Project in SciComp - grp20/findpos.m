function index=findpos(vector,val)
    if val < vector(1) || val > vector(end)
        index = -1;
    else
        for i = 2:length(vector)
            if val > vector(i-1) && val < vector(i)
                index = i;
                return;
            end
        end
    end

end
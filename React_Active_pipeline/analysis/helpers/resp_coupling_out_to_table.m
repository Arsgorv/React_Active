function T = resp_coupling_out_to_table(sessname, Out)
rows = {};
bandsLow = {'delta','theta'};
bandsMid = {'gamma','highgamma'};
for gi = 1:numel(Out.groups)
    gname = string(Out.groups{gi});

    % LOW
    if isfield(Out,'low') && numel(Out.low.group) >= gi && isfield(Out.low.group(gi),'xcorr')
        X = Out.low.group(gi).xcorr;
        % band1=delta, band2=theta
        for b = 1:2
            if b==1, W = X.band1.window; bname=bandsLow{1};
            else,    W = X.band2.window; bname=bandsLow{2};
            end
            for wi = 1:numel(W)
                s = W(wi).summary;
                rows(end+1,:) = {string(sessname), gname, "LOW", string(bname), string(W(wi).name), ...
                    s.nUsed, s.medianLag, s.iqrLag, s.meanR, s.medianR}; %#ok<AGROW>
            end
        end
    end

    % MID
    if isfield(Out,'mid') && numel(Out.mid.group) >= gi && isfield(Out.mid.group(gi),'xcorr')
        X = Out.mid.group(gi).xcorr;
        for b = 1:2
            if b==1, W = X.band1.window; bname=bandsMid{1};
            else,    W = X.band2.window; bname=bandsMid{2};
            end
            for wi = 1:numel(W)
                s = W(wi).summary;
                rows(end+1,:) = {string(sessname), gname, "MID", string(bname), string(W(wi).name), ...
                    s.nUsed, s.medianLag, s.iqrLag, s.meanR, s.medianR}; %#ok<AGROW>
            end
        end
    end
end
if isempty(rows)
    rows = cell(0,10); % must match number of VariableNames
end
T = cell2table(rows, 'VariableNames', ...
    {'session','group','spec','band','window','nUsed','medianLag_s','iqrLag_s','meanR','medianR'});
end

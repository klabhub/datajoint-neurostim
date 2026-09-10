function EEG =removeNan(EEG)
 nanData = all(isnan(EEG.data),1);
                if any(nanData(:))
                    d = diff(nanData);
                    start  =find([nanData(1) d]>0)';
                    stop = find([d<0 nanData(end)])';
                    EEG = eeg_eegrej(EEG,[start stop]);
                end
end
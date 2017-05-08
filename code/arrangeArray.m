function recChannel = arrangeArray(tdt,rec)
% This function deals with what was read from the buffer and splits it
% into nChan channels of equal length - the trick being that when half a
% buffer was read, the first value did not necessarily correspond to the
% first channel. Normally it is because we chose sizeBuf such that 
% tdt.sizeBuf/(tdt.nb_channels*2) be an integer.
% 
% Inputs: 
%   - tdt is a structure containing various parameters of our
% experiment
%   - rec is a cell containing single arrays. Each array contains a
%   sequence of doubles that were read from RZ2's buffer
% Outputs: 
%   - recChannel is a nChan*k double array where
%     - cChan is the number of channels used in the experiment
%     - k is the number of values recorded by each channel during the expt


nChan = tdt.nb_channels;
recChannel = cell(nChan,1);
indFinal = 0;

for bb=1:length(rec)
  cRec = rec{bb};
  indIniti = floor(indFinal + 1);
  indFinal = floor(indFinal + length(cRec)/nChan);
  if ~(isinteger(indFinal)&&isinteger(indIniti))
   %^% error('roblem with the indices: at least one is not an integer');
  end
  
  % Identify 
  for cc=1:nChan
    recChannel{cc}(indIniti:indFinal) = cRec(cc:nChan:end); % make sure the final cycle is also a multiple of nChan
  end
end


% figure;
% for kk=1:1:128
%   subplot(8,16,kk);
%   r = recChannel{kk};
%   plot(r); %33 34 :128:end
%   title(sprintf('kk=%d', kk));
% end

end
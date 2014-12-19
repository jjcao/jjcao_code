function cn = channels_transfer(ct, cs)
% ct: target channels, ct(2,:) is 2nd target channel
% cr: source channels
% cn: new channels
for i = 1:size(ct, 1)
    cn(i,:) = (std(ct(i,:))/std(cs(i,:))) * (cs(i,:) - mean(cs(i,:))) + mean(ct(i,:));
end

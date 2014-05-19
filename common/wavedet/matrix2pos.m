function [s] = matrix2pos(pos)

if size(pos,2) == 10
    s.Pon = pos(:,1)';
    s.P = pos(:,2)';
    s.Poff = pos(:,3)';
    s.QRSon = pos(:,4)';
    s.qrs = pos(:,5)';
    s.QRSoff = pos(:,6)';
    s.Ton = pos(:,7)';
    s.T = pos(:,8)';
    s.Tprima = pos(:,9)';
    s.Toff = pos(:,10)';
    s.Ttipo = nan(size(pos(:,10)',1),size(pos(:,10)',2));
else
    s.Pon = pos(:,1)';
    s.P = pos(:,2)';    
    s.Poff = pos(:,4)';
    s.qrs = pos(:,5)';
    s.QRSon = pos(:,6)';        
    s.QRSoff = pos(:,11)';
    s.Ton = pos(:,12)';
    s.T = pos(:,13)';
    s.Tprima = pos(:,14)';
    s.Toff = pos(:,15)';
    s.Ttipo = nan(size(pos(:,15)',1),size(pos(:,15)',2));
end
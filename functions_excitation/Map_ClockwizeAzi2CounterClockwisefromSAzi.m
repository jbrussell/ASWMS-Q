function [ AZI_DT ] = Map_ClockwizeAzi2CounterClockwisefromSAzi( AZI )
% Convert clockwise azimuth from north to Dahlen and Tromp convention,
% measured counterclockwise from due south...

 AZI_DT = AZI;
 AZI_DT(AZI <= 180) = 180-AZI_DT(AZI <= 180);
 AZI_DT(AZI > 180) = 360- (AZI_DT(AZI > 180) - 180);
 %AZI = AZI_DT;

end


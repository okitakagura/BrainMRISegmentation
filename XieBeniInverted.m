function xieBeni = XieBeniInverted(U, uExpoent, center, data)
%XIEBENIINVERTED Implementation of the measure of cluster validation of the Xie-Beni.
%   xieBeni = XieBeniInverted(U, uExpoent, center, data)
%
%   Implementation of the following article about the cluster validation
%   measures: X.L. Xie, G. Beni, "A Validity Measure for Fuzzy Clustering,"
%   IEEE Transactions on Pattern Analysis and Machine Intelligence, pp.
%   841-847, August, 1991.
%
%       See also FCMFPDEMO, INITFCM, IRISFCM, STEPFCMFP, ITERFCMFP, DISTFCMFP, and FCMFP.

%	Silvio Filipe, 31-08-2011.
%       $Revision: 2.00 $  $Date: 2011/08/31 $

% Calculate Xie-Beni Inverted cluster validation method

dist = sum(sum((U.^uExpoent) .* (distfcm(center, data).^2)));
minimum = intmax;
for i = 1 : size(center, 1)-1
    minimumAux = min(distfcmfp(center(i, :), center(i+1:end, :)).^2);
    if(minimum > minimumAux)
        minimum = minimumAux;
    end
end
xieBeni = (size(data, 1) .* minimum) / dist;

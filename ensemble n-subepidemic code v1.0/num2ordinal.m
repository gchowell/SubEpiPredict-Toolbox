function cad1=num2ordinal(num1)
%NUM2ORDINAL Ordinal text of a positive integer: 1st, 2nd, 3rd, 4th, ..., 11th, 21st.
%   cad1 = num2ordinal(num1)
%
%   Used in the figure titles of the ranked models. The earlier version listed
%   the cases 1 to 10 and returned nothing for any other number, so the
%   plotting scripts stopped at the 11th-ranked model. For 1 to 10 the text is
%   the same as before.

if ~(isnumeric(num1) && isreal(num1) && isscalar(num1) && isfinite(num1) && num1>=1 && num1==fix(num1))
    error('SubEpiPredict:InvalidRank','num2ordinal needs one positive integer.');
end
num1=full(double(num1));

lastTwo=mod(num1,100);
if lastTwo>=11 && lastTwo<=13
    suffix='th'; % 11th, 12th, 13th, 111th
else
    switch mod(num1,10)
        case 1
            suffix='st';
        case 2
            suffix='nd';
        case 3
            suffix='rd';
        otherwise
            suffix='th';
    end
end

cad1=sprintf('%d%s',num1,suffix);
end

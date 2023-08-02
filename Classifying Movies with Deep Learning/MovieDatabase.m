function d = MovieDatabase( file )

if( nargin < 1 )
  Demo
  return
end

f = fopen(file,'r');

d.name    = {};
d.rating  = [];
d.length  = [];
d.genre   = {};
d.mPAA    = {};
t         = sprintf('\t');  % a tab character
k         = 0;

while(~feof(f))
  k             = k + 1;
  q             = fgetl(f);      % one line of the file
  j             = strfind(q,t);  % find the tabs in the line
  d.name{k}     = q(1:j(1)-1);   % the name is the first token
  d.rating(1,k)	= str2double(q(j(1)+1:j(2)-1));
  d.genre{k}    = q(j(2)+1:j(3)-1);
	d.length(1,k)	= str2double(q(j(3)+1:j(4)-1));
  d.mPAA{k}     = q(j(4)+1:end);
end % end of the file

if( max(d.rating) == 0 || isnan(d.rating(1)) )
  d.rating = randi(5,1,k);
end

if( max(d.length) == 0 || isnan(d.length(1)))
  d.length = 1.8 + 0.15*randn(1,k);
end

fclose(f);


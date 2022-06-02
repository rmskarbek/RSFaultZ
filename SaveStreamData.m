function varargout = SaveStreamData(fn,varargin)
% Functions to save a stream of data to a file. The data should be
% column-oriented, so that one is building a wide matrix by appending colums.
% 
% ssd = Init(filename,[varargin])
%   Save a streaming matrix having nr rows to a designated file.
%     varargin (parameter-value pairs):
%     'nr': Number of rows in data. If not provided, it will be deduced in
%   Write. Note that this means the data matrix passed to write must be
%   sized correctly.
%     'bfr_ncalls': Store in memory for ncalls to Write before writing to
%   disk.
%
% ssd = Write(ssd,data,[force_write])
%   data is a matrix having ncols columns. If force_write is provided and is
% true, then force a write despite buffering.
%
% Finalize(ssd)
%   Write any buffered data and close the ssd struct.
%
% data = Read(filename,[cidxs],[ridxs],[isStride])
%   Read in the contents of a SSD file. If the optional argument cidxs is
%   provided, then data contains the columns listed in cidxs in ascending
%   order. cidxs can contain column numbers that exceed the number of rows in
%   the file. [ridxs] does the equivalent for rows. Set cidxs=[] if you want
%   to specifiy ridxs but not cidxs. If isStride=1 is given as the fourth
%   argument, cidxs is interpreted as a stride and should be a scalar. Set
%   ridxs=[] if you want to specifiy isStride but not ridxs.
%
% See the function 'Test' below for an example.

  [varargout{1:nargout}] = feval(fn,varargin{:});
  
function ssd = Init(fn,varargin)
  [nr bfr_ncalls] = process_options(varargin,'nr',[],'bfr_ncalls',0);
  fid = fopen(fn,'r');
  if(false && fid ~= -1)
    fclose(fid);
    error([fn ' already exists; I don''t want to overwrite it.']);
  end
  ssd.fid = fopen(fn,'wb','ieee-le');
  ssd.nr = nr;
  fwrite(ssd.fid,nr,'int32');
  ssd.bfr.ncalls = bfr_ncalls;
  ssd.bfr.i = 1;
  ssd.bfr.data = [];
  
function ssd = Write(ssd,data,force)
  if(isempty(ssd.nr))
    ssd.nr = size(data,1);
    if(ftell(ssd.fid) == 0)
      fwrite(ssd.fid,ssd.nr,'int32');
    end
  end
  if(ssd.bfr.ncalls)
    if(nargin < 3) force = 0; end
    write = force;
    ssd.bfr.data = [ssd.bfr.data data];
    if(write)
      ssd.bfr.i = 1;
    else
      if(ssd.bfr.i == ssd.bfr.ncalls)
	write = 1;
	ssd.bfr.i = 1;
      else
	ssd.bfr.i = ssd.bfr.i + 1;
      end
    end
    if(write)
      fwrite(ssd.fid,ssd.bfr.data(:),'double');
      ssd.bfr.data = [];
    end
  else
    fwrite(ssd.fid,data(:),'double');
  end
  
function Finalize(ssd)
  if(ssd.bfr.ncalls) ssd = Write(ssd,[],1); end
  fclose(ssd.fid);

function nr = NbrRows(fn)
  fid = Fopen(fn);
  nr = fread(fid,1,'int32');
  fclose(fid);
  
function fid = Fopen(fn)
  fid = fopen(fn,'rb','ieee-le');
  if(fid == -1)
    fid = fopen([fn '.dat'],'rb','ieee-le');
  end
  if(fid == -1) error(sprintf('Can''t read %s',fn)); end
  
function A = Read(fn,cidxs,ridxs,doStride)
  if(nargin < 4) doStride = 0; end
  if(nargin < 3) ridxs = []; end
  if(nargin < 2) cidxs = []; end
  doStride = doStride || isempty(cidxs);
  fid = Fopen(fn);
  A = ReadFid(fid,cidxs,ridxs,doStride);
  fclose(fid);
  
function A = Copy(src,dst,cidxs)
  fsrc = Fopen(fsrc);
  ssd = Init(dst);
  for(i = 1:length(cidxs))
    d = ReadFid(fsrc,cidxs(i));
    ssd = Write(ssd,d);
  end
  fclose(fsrc);
  Finalize(ssd);
  
function A = ReadFid(fid,cidxs,ridxs,doStride)
  if(isempty(cidxs) && isempty(ridxs))
    % Read the whole file
    nr = fread(fid,1,  'int32');
    A  = fread(fid,inf,'double');
    nc = length(A)/nr;
    A = reshape(A,nr,nc);
  elseif(doStride)
    % Read cols with a stride and possibly selected rows
    nr = fread(fid,1,'int32')';
    if(isempty(ridxs)) ridxs = 1:nr; end
    stride = cidxs;
    if(length(stride) > 1) error('stride is a scalar'); end
    if(isempty(stride)) stride = 1; end
    chunk = 1000;
    A = zeros(length(ridxs),chunk);
    cnt = 1;
    while(true)
      data = fread(fid,nr,'double');
      % Check for EOF
      if(length(data) ~= nr)
	% Resize A if necessary
	A = A(:,1:cnt-1);
	break;
      end
      % Alloc more for A?
      if(size(A,2) == cnt)
	A = [A A];
      end
      % Save this column
      A(:,cnt) = data(ridxs);
      cnt = cnt + 1;
      if(stride > 1)
	% Skip ahead
	stat = fseek(fid,8*(stride-1)*nr,'cof');
	if(stat < 0)
	  A = A(:,1:cnt-1);
	  break;
	end
      end
    end
  else
    % Read selected cols and possibly rows
    nr = fread(fid,1,'int32')';
    if(isempty(ridxs)) ridxs = 1:nr; end
    ni = length(cidxs);
    A = zeros(length(ridxs),ni);
    d = 0;
    for(i = 1:length(cidxs))
      % Nbr of cols between next desired one and the one we're reading now
      d = cidxs(i) - d - 1;
      if(d > 0)
	% Skip ahead. Unfortunately, I have to hardcode 8 bytes, as fseek
	% does not except a precision argument as do fread and fwrite.
	stat = fseek(fid,8*d*nr,'cof');
	if(stat < 0)
	  A = A(:,1:i-1);
	  break;
	end 
      end
      data = fread(fid,nr,'double');
      if(length(data) ~= nr)
	% Resize A if necessary
	A = A(:,1:i-1);
	break;
      end
      A(:,i) = data(ridxs);
      d = cidxs(i);
    end
  end
  
function Test
  fn = 'ssdtest.dat';
  ssd = SaveStreamData('Init',fn);
  A0 = (1:10)'*(1:20);
  SaveStreamData('Write',ssd,A0);
  SaveStreamData('Finalize',ssd);
  A1 = SaveStreamData('Read',fn);
  A2 = SaveStreamData('Read',fn,2:2:100);
  A0,A1,A2

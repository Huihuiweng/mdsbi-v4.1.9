function data = loadmesh(pb,name,st,tw)
  
% load data
  
  eval(['ms = pb.',name,';'])
  
  if isempty(ms)
    disp('field not defined')
    return
  end
    
  if strcmp(ms.name,'front')==0
  
    % read base time vector
    
    [f m] = fopen([pb.f00,'_',ms.name,'.t'],'r',pb.endian);
    if f==-1
      if ~isempty(m);disp(m);return,end
    else
      eval(['t = fread(f,',pb.pre,');']);
      n = length(t);
      fclose(f);
    end

  else
    
    n = 1;
    t = 0;
    
  end
    
  disp(['Total time steps = ',num2str(n)])
  
  % assign position vectors
  
  data = ms;
  
  data.x = linspace(ms.xmin,ms.xmax,ms.nx);
  data.y = linspace(ms.ymin,ms.ymax,ms.ny);
  data.z = linspace(ms.zmin,ms.zmax,ms.nz);

  % strided and windowed times to load
  
  if nargin==2; 
    st = 1; % factor multiplying stride in time
    tw = [min(t) max(t)]; % time window
  elseif nargin==3
    tw = [min(t) max(t)];
  end
  
  tw(1) = max(ms.tmin,tw(1));
  nmin = find(t>=tw(1),1,'first');

  tw(2) = min(ms.tmax,tw(2));
  nmax = find(t<=tw(2),1,'last');
  
  data.t = t(nmin:st:nmax);
  
  data.nt = length(data.t);
  data.tmin = min(data.t);
  data.tmax = max(data.t);
  
  % set dimensions
  
  dim = 3;
  if ms.nx==1;dim=dim-1;end
  if ms.ny==1;dim=dim-1;end
  if ms.nz==1;dim=dim-1;end

  % open file
  
  [f m] = fopen([pb.f00,'_',ms.name,'.dat'],'r',pb.endian);
  if ~isempty(m);disp(m);end
  if f==-1;return;end
  
  % indices of time steps to be loaded

  if length(t)==data.nt
    it = 1:length(t);
  else
    it = zeros(1,data.nt);
    for n=1:data.nt
      it(n) = find(data.t(n)==t);
    end
  end
      
  % read data
  
  switch dim
    
   case 0
    
    eval(['data.',ms.field,' = zeros(1,data.nt);']);
    readdata = ['data.',ms.field,'(n) = fread(f,1,',pb.pre,');'];
    datasize = 1;
    
   case 1
    
    if ms.nx~=1

      eval(['data.',ms.field,' = zeros(ms.nx,data.nt);']);
      datasize = ms.nx;
      readdata = ['data.',ms.field,'(:,n) = fread(f,ms.nx,',pb.pre,');'];
    
    elseif ms.ny~=1
      
      eval(['data.',ms.field,' = zeros(ms.ny,data.nt);']);
      readdata = ['data.',ms.field,'(:,n) = fread(f,ms.ny,',pb.pre,');'];
      datasize = ms.ny;
    
    elseif ms.nz~=1

      eval(['data.',ms.field,' = zeros(ms.nz,data.nt);']);
      readdata = ['data.',ms.field,'(:,n) = fread(f,ms.nz,',pb.pre,');'];
      datasize = ms.nz;
    
    end
   
   case 2

    if ms.nx==1

      eval(['data.',ms.field,' = zeros(ms.ny,ms.nz,data.nt);']);
      readdata = ['data.',ms.field,'(:,:,n) = fread(f,[ms.ny ms.nz],',pb.pre,');'];
      datasize = ms.ny*ms.nz;
    
    elseif ms.ny==1

      eval(['data.',ms.field,' = zeros(ms.nx,ms.nz,data.nt);']);
      readdata = ['data.',ms.field,'(:,:,n) = fread(f,[ms.nx ms.nz],',pb.pre,');'];
      datasize = ms.nx*ms.nz;
    
    elseif ms.nz==1

      eval(['data.',ms.field,' = zeros(ms.nx,ms.ny,data.nt);']);
      readdata = ['data.',ms.field,'(:,:,n) = fread(f,[ms.nx ms.ny],',pb.pre,');'];
      datasize = ms.nx*ms.ny;
    
    end
   
   case 3

    eval(['data.',ms.field,' = zeros(ms.nx,ms.ny,ms.nz,data.nt);']);
    readdata = ['for iz=1:ms.nz; data.',ms.field,'(:,:,iz,n) = fread(f,[ms.nx ms.ny],',pb.pre,'); end'];
    datasize = ms.nx*ms.ny*ms.nz;
  
  end
  
  disp(['loading ',num2str(data.nt),' time steps'])
  it_old = 0;
  for n=1:data.nt
    nskip = it(n)-it_old-1;
    offset = 4*datasize*nskip;
    fseek(f,offset,0);
    eval(readdata)
    it_old = it(n);
  end
      
  fclose(f);
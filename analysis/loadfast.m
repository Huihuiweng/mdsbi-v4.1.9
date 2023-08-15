function data = loadfast(pb,name,nt,st)
  
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
      if nargin>=3
        eval(['t = fread(f,nt,',pb.pre,');']);
      else
        eval(['t = fread(f,',pb.pre,');']);
        nt = length(t);
      end
      fclose(f);
    end

  else
    
    nt = 1;
    data.t = 0;
    
  end

  disp(['Total time steps = ',num2str(length(t))])
  
  % assign position vectors
  
  data = ms;
  
  data.x = linspace(ms.xmin,ms.xmax,ms.nx);
  data.y = linspace(ms.ymin,ms.ymax,ms.ny);
  data.z = linspace(ms.zmin,ms.zmax,ms.nz);

  data.t = t;
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
  
  % read data
  
  switch dim
    
   case 0
    
    readdata = ['data.',ms.field,' = fread(f,nt,',pb.pre,');'];
    
   case 1
    
    if ms.nx~=1

      eval(['data.',ms.field,' = zeros(ms.nx,data.nt);']);
      readdata = ['data.',ms.field,' = fread(f,[ms.nx nt],',pb.pre,');'];
    
    elseif ms.ny~=1
      
      eval(['data.',ms.field,' = zeros(ms.ny,data.nt);']);
      readdata = ['data.',ms.field,' = fread(f,[ms.ny nt],',pb.pre,');'];
    
    elseif ms.nz~=1

      eval(['data.',ms.field,' = zeros(ms.nz,data.nt);']);
      readdata = ['data.',ms.field,' = fread(f,[ms.nz nt],',pb.pre,');'];
    
    end
   
   case 2

    if ms.nx==1

      eval(['data.',ms.field,' = zeros(ms.ny,ms.nz,data.nt);']);
      readdata = ['data.',ms.field,'(:,:,n) = fread(f,[ms.ny ms.nz],',pb.pre,');'];
    
    elseif ms.ny==1

      eval(['data.',ms.field,' = zeros(ms.nx,ms.nz,data.nt);']);
      readdata = ['data.',ms.field,'(:,:,n) = fread(f,[ms.nx ms.nz],',pb.pre,');'];
    
    elseif ms.nz==1

      eval(['data.',ms.field,' = zeros(ms.nx,ms.ny,data.nt);']);
      readdata = ['data.',ms.field,'(:,:,n) = fread(f,[ms.nx ms.ny],',pb.pre,');'];
    
    end
   
   case 3

    eval(['data.',ms.field,' = zeros(ms.nx,ms.ny,ms.nz,data.nt);']);
    readdata = ['data.',ms.field,'(:,:,k,n) = fread(f,[ms.nx ms.ny],',pb.pre,');'];
    
  end

  disp(['loading ',num2str(data.nt),' time steps'])

  switch dim
   case 2
    for n=1:data.nt
      eval(readdata)
    end
   case 3
    for n=1:data.nt
      for k=1:ms.nz
	eval(readdata)
      end
    end
   otherwise
    eval(readdata) 
  end
  
  fclose(f);
program imdraw1;
uses crt,graph,dos;
{-------------------------------
 Gabor and Tamas Szabados
  1996-97
 -------------------------------}
const
  n=600000;
  m=60000;
  ncoord=4;
type vectn=array[1..n] of integer;
     vectm=array[1..m] of integer;
     vectc=array[1..ncoord] of integer;
     vectmc=array[1..m,1..ncoord] of integer;
var
  szin,valasz,lehet:byte;
  errcode,i,j,tpill,bgcolor:integer;
  db,lep,kesl,tkezdo,tbef,t1k,t1b:integer;
  db2,db3,tlep,kesl2,ugr,tkezdo2,tbef2,ylep,telso,tutso:integer;
  nwtyp,nrtyp:integer;
  c,megjel,arnyek:char;
  tim:string[5];
  inp:text;
  xw,yw,xr,yr:^vectc;
  x,y,r,t1,t2,xt,yt,rt,tht1,tht2,ind,indt:^vectn;
  time,bcells,marrow,abody,thelper,ilns:^vectm;
  self,nonself:^vectmc;
label
  hiba;

{---------------------------------------------------------------------------}
procedure kezdoertek;
begin
  for i:=1 to n do
    begin
      x^[i]:=0;
      y^[i]:=0;
      r^[i]:=0;
      t1^[i]:=0;
      t2^[i]:=0;
      xt^[i]:=0;
      yt^[i]:=0;
      rt^[i]:=0;
      tht1^[i]:=0;
      tht2^[i]:=0;
      ind^[i]:=0;
      indt^[i]:=0;
      if i<=m then
        begin
          time^[i]:=0;
          bcells^[i]:=0;
          marrow^[i]:=0;
          abody^[i]:=0;
          thelper^[i]:=0;
          ilns^[i]:=0;
          for j:=1 to ncoord do
            begin
              self^[i,j]:=0;
              nonself^[i,j]:=0;
            end{j};
        end;
    end;
  szin:=textattr;
end;

{---------------------------------------------------------------------------}
procedure adatbe1;
var
  dirinfo:searchrec;
  fajl:array[1..50] of string[30];
begin
  i:=0;
  j:=0;
  clrscr;
  findfirst('o*.*',archive,dirinfo);
  while doserror=0 do
    begin
      inc(i);
      writeln(i,'. ',dirinfo.name);
      fajl[i]:=''+dirinfo.name;
      findnext(dirinfo);
    end;
  writeln;
  writeln('Choose an output file by typing in its number:');
  repeat
    gotoxy(57,wherey-1); clreol;
    readln(valasz);
  until (valasz>=1) and (valasz<=i);
  inc(j);
  while (j<=i) and (j<=valasz) do
    begin
      if valasz=j then assign(inp,fajl[j]);
      inc(j);
    end;
end;

{---------------------------------------------------------------------------}
procedure valaszt;
begin
  clrscr;
  normvideo;
  gotoxy(1,20);
  writeln('1. Time diagram');
  writeln('2. Shape space process');
  {writeln('3. Time diagram and shape space process together');}
  writeln('4. Shape space figure');
  writeln('5. Peptid space process');
  writeln;
  writeln('Choice:');
  repeat
    gotoxy(12,wherey-1); clreol;
    readln(valasz);
  until (valasz>=1) and (valasz<=5);
end;

{---------------------------------------------------------------------------}
procedure menu;
begin
  writeln('1. actual file again');
  writeln('2. new output file');
  writeln('3. exit');
  textattr:=szin;
  writeln;
  writeln('Choice:');
  repeat
    gotoxy(12,wherey-1); clreol;
    readln(lehet);
  until (lehet=1) or (lehet=2) or (lehet=3);
end;

{---------------------------------------------------------------------------}
procedure beolvas1;
var
  i:integer;
  be1,be2:string[10];
  be3:string[7];
  be4,be5:string[4];
  be:string[6];
begin
  db2:=0;
  be:='';
  reset(inp);
  while be<>' kth0=' do
    readln(inp,be);
  // readln(inp);
  readln(inp,be1,nwtyp,be2,nrtyp);

  for i:=1 to nwtyp do
    readln(inp,be3,be4,xw^[i],be5,yw^[i]);
  for i:=1 to nrtyp do
    readln(inp,be3,be4,xr^[i],be5,yr^[i]);

end;

{---------------------------------------------------------------------------}
procedure beolvas2;
var
  tbee:integer;
  be:string[3];
  be1:string[15];
  be2:string[2];
  t10,t20,t11,t22:real;
  sel12:string[12];
  sel5:string[5];
  sel6:string[6];
  sel38:string[38];
label
  sok;
begin
  db:=0;
  db2:=0;
  db3:=0;
  be:='';
  be1:='';
  reset(inp);
  while be<>'The' do
    readln(inp,be);
  //readln(inp,be);

  readln(inp,be);
  readln(inp,be);
  read(inp,be);
  while be<>'' do
    begin
      if be = 'bla' then
        begin
           inc(db);
          if db > n then
            begin
              errcode:=1;
              goto sok;
            end;
          readln(inp,sel12,x^[db],y^[db],r^[db],tbee,tbee,tbee,t10,t20);
          t1^[db]:=round(t10);
          t2^[db]:=abs(round(t20));
          ind^[db]:=db;
        end

      else
        begin
        if (be<>'The') then
         begin
          inc(db2);
          if db2 > m then
            begin
              errcode:=2;
              goto sok;
            end;
          read(inp,time^[db2],sel6);
          for j:=1 to 1 do
            read(inp,self^[db2,j]);

         // !!!!!!!
          read(inp,sel6);
          for j:=1 to 1 do
            read(inp,nonself^[db2,j]);
            read(inp,sel5,bcells^[db2],sel6,abody^[db2],sel6,thelper^[db2],
            sel6,ilns^[db2]);
            readln(inp,sel5,marrow^[db2]);
         // !!!!!!!
         // readln(inp,sel6,bcells^[db2],sel6,marrow^[db2]);
         end
         else if be = 'The' then
                begin
                  inc(db3);
                  if db > n then
                    begin
                    errcode:=1;
                    goto sok;
                  end;
                  readln(inp,sel6,xt^[db3],yt^[db3],rt^[db3],t11,t22);
                  tht1^[db3]:=round(t11);
                  tht2^[db3]:=abs(round(t22));
                  indt^[db3]:=db3;

                end
              else readln(inp,be);


        end;
      read(inp,be);
    end;

  while be1<>' The list of bl' do
    readln(inp,be1);
  readln(inp,be1);
  readln(inp,be1);
  read(inp,be1);
  while be1<>'' do
    begin
      inc(db);
      readln(inp,x^[db],y^[db],r^[db],tbee,tbee,tbee,t10,t20);
      t1^[db]:=round(t10);
      t2^[db]:=abs(round(t20));
      ind^[db]:=db;
      read(inp,be1);
    end;

  while be2<>' T' do //
    readln(inp,be2);

    readln(inp,be2);

    while be2<>' T' do //
    readln(inp,be2);

  readln(inp,be2);
  readln(inp,be2);
  read(inp,be2);
  while be2<>'' do
    begin
      inc(db3);
      readln(inp,xt^[db3],yt^[db3],rt^[db3],t11,t22);
      tht1^[db3]:=round(t11);
      tht2^[db3]:=abs(round(t22));
      indt^[db3]:=db3;
      read(inp,be2);
    end;


  t1k:=t1^[ind^[1]];
  t1b:=t1^[ind^[db-1]];
  telso:=time^[1];
  tutso:=time^[db2];
  sok:;
end;

{---------------------------------------------------------------------------}
procedure ellenor;

 begin
   ugr:=0;
   while x^[ugr]<> 0 do
    begin
     writeln(ugr,'  ',x^[ugr],'  ',y^[ugr],'   ',r^[ugr],'  ',ind^[ugr]);
     inc(ugr);
    end;
    delay(5000);

   ugr:=0;
    while xt^[ugr]<> 0 do
     begin
      writeln(tht1^[ugr],'   ',xt^[ugr],'  ',yt^[ugr],'   ',rt^[ugr],'   ',
      indt^[ugr]);
      inc(ugr);
     end;
    delay(10000);
     clreol;
 end;
{---------------------------------------------------------------------------}
procedure adatbe2;
begin
  lep:=0;
  textattr:=szin;
  clrscr;
  writeln('First time instant between ',telso,' and ',tutso,' : ');
  writeln;
  repeat
    gotoxy(1,wherey-1); clreol;
    readln(tkezdo);
  until (tkezdo>=telso) and (tkezdo<=tutso);
  writeln;
  writeln('Last time instant between ',tkezdo,' and ',tutso,' : ');
  writeln;
  repeat
    gotoxy(1,wherey-1); clreol;
    readln(tbef);
  until (tbef>=tkezdo) and (tbef<=tutso);
  writeln;
  writeln('Plotting with (p)oints, (c)ircles, or (s)quares ?');
  writeln;
  repeat
    gotoxy(1,wherey-1); clreol;
    readln(megjel);
  until (megjel='p') or (megjel='c') or (megjel='s');
  writeln;
  writeln('Time step:');
  writeln;
  repeat
    gotoxy(1,wherey-1); clreol;
    readln(lep);
  until (lep>0) and (lep<=tbef-tkezdo);
  writeln;
  writeln('Delay:');
  writeln;
  repeat
    gotoxy(1,wherey-1); clreol;
    readln(kesl);
  until kesl>=0;
end;

{---------------------------------------------------------------------------}
procedure adatbe3;
begin
  textattr:=szin;
  clrscr;
  writeln('The first time instant between 0 and ',tutso,': ');
  writeln;
  repeat
    gotoxy(1,wherey-1); clreol;
    readln(tkezdo2);
  until (tkezdo2>=0) and (tkezdo2<=tutso);
  writeln;
  writeln('The last time instant between ',tkezdo2,' and ',tutso,': ');
  writeln;
  repeat
    gotoxy(1,wherey-1); clreol;
    readln(tbef2);
  until (tbef2>=tkezdo2) and (tbef2<=tutso);
  writeln;
  writeln('Time step:');
  writeln;
  repeat
    gotoxy(1,wherey-1); clreol;
     readln(tlep);
  until tlep>=1;
  writeln;
  writeln('Vertical unit:');
  writeln;
  repeat
    gotoxy(1,wherey-1); clreol;
    readln(ylep);
  until ylep>=1;
  writeln;
  writeln('Delay:');
  writeln;
  repeat
    gotoxy(1,wherey-1); clreol;
    readln(kesl2);
  until kesl2>=0;

end;

{---------------------------------------------------------------------------}
procedure adatbe4;
begin
  lep:=0;
  textattr:=szin;
  clrscr;
  writeln('Time instant between ',telso,' and ',tutso,' : ');
  writeln;
  repeat
    gotoxy(1,wherey-1); clreol;
    readln(tkezdo);
  until (tkezdo>=telso) and (tkezdo<=tutso);
  writeln;
  writeln('Plotting with (p)oints, (c)ircles, or (s)quares ?');
  writeln;
  repeat
    gotoxy(1,wherey-1); clreol;
    readln(megjel);
  until (megjel='p') or (megjel='c') or (megjel='s');
  writeln;
  writeln('Shade the covered region? (y/n)');
  writeln;
  repeat
    gotoxy(1,wherey-1); clreol;
    readln(arnyek);
  until (arnyek='y') or (arnyek='n');
  writeln;
end;


{---------------------------------------------------------------------------}

procedure rendezt;
var
  mint,segedt,tht1m,tht1a:integer;
begin
  clrscr;
  textcolor(red+blink);
  writeln('Working...');
  for i:=1 to db3-1 do
    begin
      mint:=i;
      tht1m:=tht1^[indt^[mint]];
      for j:=i+1 to db3 do
        begin
          tht1a:=tht1^[indt^[j]];
          if tht1a<tht1m then
            begin
              mint:=j;
              tht1m:=tht1a;
            end;
        end;
      if mint<>i then
        begin
          segedt:=indt^[i];
          indt^[i]:=indt^[mint];
          indt^[mint]:=segedt;
        end;
    end;
end;

{---------------------------------------------------------------------------}
procedure rendez;
var
  min,seged,t1m,t1a:integer;
begin
  clrscr;
  textcolor(red+blink);
  writeln('Working...');
  for i:=1 to db-1 do
    begin
      min:=i;
      t1m:=t1^[ind^[min]];
      for j:=i+1 to db do
        begin
          t1a:=t1^[ind^[j]];
          if t1a<t1m then
            begin
              min:=j;
              t1m:=t1a;
            end;
        end;
      if min<>i then
        begin
          seged:=ind^[i];
          ind^[i]:=ind^[min];
          ind^[min]:=seged;
        end;
    end;
end;

{---------------------------------------------------------------------------}
procedure grafikon;
var
  gd,gm:smallint;
  ErrCode:integer;
begin
  gd:=detect;
  initgraph(gd,gm,'c:\bin\tp\bgi');

  ErrCode := GraphResult;
  if ErrCode = grOk then
   cleardevice
  else
   WriteLn('Graphics error:',
           GraphErrorMsg(ErrCode));

  if bgcolor=1 then
   begin
    cleardevice;
    setbkcolor(15);
    cleardevice;
   end;

end;

{---------------------------------------------------------------------------}
procedure keret;
begin
  if bgcolor = 1 then
   setcolor(16)
  else
   setcolor(white);

  rectangle(70,0,569,479);
  setlinestyle(0,0,3);
  line(70,240,569,240);
  line(70,0,70,479);
  setlinestyle(0,0,1);
  {line(320,0,320,479);}
end;

{---------------------------------------------------------------------------}
procedure negyzet(szul:integer);
var
  xx,yy,rr:integer;
  square : array[1..4] of PointType;
begin
  xx:=70+x^[ind^[i]] div 2;
  yy:=240+y^[ind^[i]] div 2;
  rr:=r^[ind^[i]] div 2;
  if szul=1 then
    begin
      if arnyek = 'y' then
        begin
          square[1].x:=xx-rr;
          square[1].y:=yy-rr;
          square[2].x:=xx+rr;
          square[2].y:=yy-rr;
          square[3].x:=xx+rr;
          square[3].y:=yy+rr;
          square[4].x:=xx-rr;
          square[4].y:=yy+rr;
          setfillstyle(1,1);
          FillPoly(SizeOf(square) div
            SizeOf(PointType), square);
        end
      else
        begin
          setcolor(green);
          rectangle(xx-rr,yy-rr,xx+rr,yy+rr);
        end;
    end
  else
    begin
      if bgcolor=1 then
       setcolor(white)
      else
       setcolor(16);

      rectangle(xx-rr,yy-rr,xx+rr,yy+rr);
      keret;
    end;
end;

{---------------------------------------------------------------------------}
procedure kor(szul:integer);
begin
  if szul=1 then
    begin
      if arnyek = 'y' then
        begin
          setfillstyle(1,1);
          FillEllipse(70+x^[ind^[i]] div 2,240+y^[ind^[i]] div 2,r^[ind^[i]] div 2,
            r^[ind^[i]] div 2);
        end
      else
        begin
          setcolor(green);
          circle(70+x^[ind^[i]] div 2,240+y^[ind^[i]] div 2,r^[ind^[i]] div 2);
        end;
    end
  else
    begin
      if bgcolor=1 then
       setcolor(white)
      else
       setcolor(16);

      circle(70+x^[ind^[i]] div 2,240+y^[ind^[i]] div 2,r^[ind^[i]] div 2);
      keret;
    end;
end;

{---------------------------------------------------------------------------}
procedure pont(szul:integer);
begin
  if szul=1 then
    begin
      if arnyek='y' then
        begin
          setcolor(blue);
          circle(70+x^[ind^[i]] div 2,240+y^[ind^[i]] div 2,1);
          setcolor(green);
          putpixel(70+x^[ind^[i]] div 2,240+y^[ind^[i]] div 2,2);
        end
      else
        begin
          setcolor(green);
          circle(70+x^[ind^[i]] div 2,240+y^[ind^[i]] div 2,1);
          putpixel(70+x^[ind^[i]] div 2,240+y^[ind^[i]] div 2,2);
        end;
    end
  else
    begin
      if bgcolor=1 then
       setcolor(white)
      else
       setcolor(16);

      circle(70+x^[ind^[i]] div 2,240+y^[ind^[i]] div 2,1);
      putpixel(70+x^[ind^[i]] div 2,240+y^[ind^[i]] div 2,0);
      keret;
    end;
end;

{---------------------------------------------------------------------------}
procedure vizsgal(szul:integer);
begin
  if megjel='p' then pont(szul)
    else if megjel='c' then kor(szul)
           else if megjel='s' then negyzet(szul);
end;

{---------------------------------------------------------------------------}
procedure thnegyzet(szul:integer);
var
  xx,yy,rr:integer;
  square : array[1..4] of PointType;
begin
  xx:=70+xt^[indt^[i]] div 2;
  yy:=240+yt^[indt^[i]] div 2;
  rr:=abs(rt^[indt^[i]]) div 2;
  if szul=1 then
    begin
      if arnyek = 'y' then
        begin
          square[1].x:=xx-rr;
          square[1].y:=yy-rr;
          square[2].x:=xx+rr;
          square[2].y:=yy-rr;
          square[3].x:=xx+rr;
          square[3].y:=yy+rr;
          square[4].x:=xx-rr;
          square[4].y:=yy+rr;
          setfillstyle(1,1);
          FillPoly(SizeOf(square) div
            SizeOf(PointType), square);
        end
      else
        begin
          setcolor(green);
          rectangle(xx-rr,yy-rr,xx+rr,yy+rr);
        end;
    end
  else
    begin
      if bgcolor=1 then
       setcolor(white)
      else
       setcolor(16);

      rectangle(xx-rr,yy-rr,xx+rr,yy+rr);
      keret;
    end;
end;

{---------------------------------------------------------------------------}
procedure thkor(szul:integer);
begin
  if szul=1 then
    begin
      if arnyek = 'y' then
        begin
          setfillstyle(1,1);
          FillEllipse(70+xt^[indt^[i]] div 2,240+yt^[indt^[i]] div 2,abs(rt^[indt^[i]]) div 2,
            abs(rt^[indt^[i]]) div 2);
        end
      else
        begin
          setcolor(green);
          circle(70+xt^[indt^[i]] div 2,240+yt^[indt^[i]] div 2,abs(rt^[indt^[i]]) div 2);
        end;
    end
  else
    begin
      if bgcolor=1 then
       setcolor(white)
      else
       setcolor(16);

      circle(70+xt^[indt^[i]] div 2,240+yt^[indt^[i]] div 2,abs(rt^[indt^[i]]) div 2);
      keret;
    end;
end;

{---------------------------------------------------------------------------}
procedure thpont(szul:integer);
begin
  if szul=1 then
    begin
      if arnyek='y' then
        begin
          setcolor(blue);
          circle(70+xt^[indt^[i]] div 2,240+yt^[indt^[i]] div 2,1);
          setcolor(green);
          putpixel(70+xt^[indt^[i]] div 2,240+yt^[indt^[i]] div 2,2);
        end
      else
        begin
          setcolor(green);
          circle(70+xt^[indt^[i]] div 2,240+yt^[indt^[i]] div 2,1);
          putpixel(70+xt^[indt^[i]] div 2,240+yt^[indt^[i]] div 2,2);
        end;
    end
  else
    begin
      if bgcolor=1 then
       setcolor(white)
      else
       setcolor(16);

      circle(70+xt^[indt^[i]] div 2,240+yt^[indt^[i]] div 2,1);
      putpixel(70+xt^[indt^[i]] div 2,240+yt^[indt^[i]] div 2,0);
      keret;
    end;
end;

{---------------------------------------------------------------------------}
procedure vizsgalth(szul:integer);
begin
  if megjel='p' then thpont(szul)
    else if megjel='c' then thkor(szul)
           else if megjel='s' then thnegyzet(szul);
end;

{---------------------------------------------------------------------------}
procedure ciklus;
var
  ii:integer;
begin
  arnyek:='n';
  tpill:=tkezdo;
  while tpill<=tbef do
    begin
      i:=1;
      while (i<=db) and (t1^[ind^[i]]<tpill) do
        begin
          if (t2^[ind^[i]]>=tpill-lep) and (t2^[ind^[i]] < tpill) then
            vizsgal(2);
          inc(i);
        end;
      i:=1;
      while (i<=db) and (t1^[ind^[i]]<=tpill) do
        begin
          if t2^[ind^[i]]>=tpill then
            vizsgal(1);
          inc(i);
        end;
      setcolor(yellow);
      for ii:=1 to nwtyp do
        outtextxy(66+xw^[ii] div 2,236+yw^[ii] div 2,'s');
      setcolor(red);
      for ii:=1 to nrtyp do
        outtextxy(66+xr^[ii] div 2,236+yr^[ii] div 2,'n');
      str(tpill:5,tim);
      setfillstyle(1,0);

      {bar(0,0,67,17);}

      if bgcolor=1 then
       setcolor(16)
      else
       setcolor(white);

      bar(0,0,69,479);
      bar(570,0,639,479);
      settextstyle(3,0,2);
      outtextxy(50,0,'y');
      outtextxy(575,224,'x');
      keret;
      settextstyle(0,0,1);
      outtextxy(0,50,concat('t= ',tim));
      delay(kesl);
      inc(tpill,lep);
      if keypressed then
        begin
          c:=readkey;
          c:=readkey;
        end;
    end;
    c:=readkey;
    closegraph;
end;


{---------------------------------------------------------------------------}
procedure ciklust;
var
  ii:integer;
begin
  arnyek:='n';
  tpill:=tkezdo;
  while tpill<=tbef do
    begin
      i:=1;
      while (i<=db3) and (tht1^[indt^[i]]<tpill) do
        begin
          if (tht2^[indt^[i]]>=tpill-lep) and (tht2^[indt^[i]] < tpill) then
            vizsgalth(2);
          inc(i);
        end;
      i:=1;
      while (i<=db3) and (tht1^[indt^[i]]<=tpill) do
        begin
          if tht2^[indt^[i]]>=tpill then
            vizsgalth(1);
          inc(i);
        end;
      setcolor(yellow);
      for ii:=1 to nwtyp do
        outtextxy(66+xw^[ii] div 2,236+yw^[ii] div 2,'s');
      setcolor(red);
      for ii:=1 to nrtyp do
        outtextxy(66+xr^[ii] div 2,236+yr^[ii] div 2,'n');
      str(tpill:5,tim);
      setfillstyle(1,0);

      {bar(0,0,67,17);}

      if bgcolor=1 then
       setcolor(16)
      else
       setcolor(white);

      bar(0,0,69,479);
      bar(570,0,639,479);
      settextstyle(3,0,2);
      outtextxy(50,0,'y');
      outtextxy(575,224,'x');
      keret;
      settextstyle(0,0,1);
      outtextxy(0,50,concat('t= ',tim));
      delay(kesl);
      inc(tpill,lep);
      if keypressed then
        begin
          c:=readkey;
          c:=readkey;
        end;
    end;
    c:=readkey;
    closegraph;
end;

{---------------------------------------------------------------------------}
procedure ciklus4;
var
 ii:integer;
begin
  if arnyek = 'y' then
    begin
      i:=1;
      while (i<=db) and (t1^[ind^[i]]<=tkezdo) do
        begin
          if t2^[ind^[i]]>=tkezdo then
          vizsgal(1);
          inc(i);
        end;
      arnyek:='n';
    end;
  i:=1;
  while (i<=db) and (t1^[ind^[i]]<=tkezdo) do
    begin
      if t2^[ind^[i]]>=tkezdo then
        vizsgal(1);
      inc(i);
    end;
  setcolor(yellow);
  for ii:=1 to nwtyp do
    outtextxy(66+xw^[ii] div 2,236+yw^[ii] div 2,'s');
  setcolor(red);
  for ii:=1 to nrtyp do
  outtextxy(66+xr^[ii] div 2,236+yr^[ii] div 2,'n');
  setfillstyle(1,0);
  bar(0,0,69,479);
  bar(570,0,639,479);
  keret;
  settextstyle(3,0,2);
  outtextxy(50,0,'y');
  outtextxy(575,224,'x');
  c:=readkey;
  closegraph;
end;

{---------------------------------------------------------------------------}
procedure keret2;
type
  nevstr=string[10];
var
  ii,alap:integer;
  szov:string[9];
procedure jel(szin:integer);
begin
  setcolor(szin); outtextxy(0,alap,'Û');
  inc(alap,20);
end;
procedure cimke(szoveg:nevstr);
begin
  outtextxy(10,alap,szoveg);
  inc(alap,20);
end;
procedure cimke2(koord:integer;szoveg:integer);
begin
  str(szoveg:9,szov);
  outtextxy(koord,470,szov);
end;
procedure cimke3(koord:integer;szoveg:integer);
begin
  str(szoveg:9,szov);
  outtextxy(0,koord,szov);
end;
begin
  alap:=80;
  jel(14);
  jel(4);
  jel(2);
  jel(6);
  jel(9);
  jel(12);
  jel(10);
  if bgcolor = 1 then
   setcolor(16)
  else
   setcolor(white);

  alap:=80;
  cimke('self');
  cimke('nonself');
  cimke('Bcells');
  cimke('marrow');
  cimke('Thcells');
  cimke('antibodys');
  cimke('i.leukins');
  rectangle(80,0,639,460);
  for i:=1 to 5 do
    begin
      line(80+i*100,457,80+i*100,463);
      line(77,460-i*100,83,460-i*100);
    end;
  cimke2(15,tkezdo2);
  for ii:=1 to 5 do
    cimke2(ii*100+25,tkezdo2+100*ii*tlep);
  outtextxy(630,470,'t');
  for ii:=1 to 4 do
    cimke3(457-ii*100,ii*ylep);
  outtextxy(56,0,'no.');
end;

{---------------------------------------------------------------------------}
procedure diagramm;
var
  i,k1,k2,i0,tt2:integer;
  id,yd,tk:real;
label
  kesz;
begin
  k1:=-1;
  id:=tlep/10;
  yd:=ylep/100;
  tk:=tkezdo2/tlep;
  for i:=1 to db2-1 do
    begin
      tt2:=time^[i];
      if (tt2 >= tkezdo2) and ((tt2 mod tlep) = 0) then
        begin
          k2:=round((tt2/tlep)-tk)+80;
          if k1 >= 0 then
            begin
              setcolor(14);
              for j:=1 to nwtyp do
                line(k1,460-round(self^[i0,j]/ yd),k2,
                     460-round(self^[i,j]/yd));
              setcolor(4);
              for j:=1 to nrtyp do
                line(k1,460-round(nonself^[i0,j]/yd),k2,
                     460-round(nonself^[i,j]/yd));
              setcolor(2);
              line(k1,460-round(bcells^[i0]/yd),k2,460-round(bcells^[i]/yd));
              setcolor(6);
              line(k1,460-round(marrow^[i0]/yd),k2,460-round(marrow^[i]/yd));
              setcolor(9);
              line(k1,460-round(thelper^[i0]/yd),k2,460-round(thelper^[i]/yd));
              setcolor(12);
              line(k1,460-round(abody^[i0]/yd),k2,460-round(abody^[i]/yd));
              setcolor(10);
              line(k1,460-round(ilns^[i0]/yd),k2,460-round(ilns^[i]/yd));
              delay(kesl2);
            end{if};
          i0:=i;
          k1:=k2;
        end{if};
      if tt2 >= tbef2 then
        goto kesz;
    end{i};
    kesz:;
    c:=readkey;
    closegraph;
end;

{***************************************************************************}
begin {foprogram kezdete}
  new(xw);
  new(yw);
  new(xr);
  new(yr);
  new(x);
  new(y);
  new(r);
  new(t1);
  new(t2);
  new(ind);
  new(indt);
  new(time);
  new(self);
  new(bcells);
  new(nonself);
  new(marrow);
  new(xt);
  new(yt);
  new(rt);
  new(tht1);
  new(tht2);
  new(abody);
  new(thelper);
  new(ilns);
  errcode:=0;
  bgcolor:=1;
  repeat
    kezdoertek;
    adatbe1;
    beolvas1;
    beolvas2;

    if errcode <> 0 then
      goto hiba;
    rendez;
    rendezt;
    repeat
     valaszt;
     if valasz=2 then
       begin
           adatbe2;
           grafikon;
           keret;
           ciklus;
           menu;
       end
       else if valasz=1 then
         begin
             adatbe3;
             grafikon;
             keret2;
             diagramm;
             menu;
         end
       else if valasz=5 then
         begin
             //ellenor;
             adatbe2;
             grafikon;
             keret;
             ciklust;
             menu;
         end
       else if valasz=4 then
         begin
             adatbe4;
             grafikon;
             keret;
             ciklus4;
             menu;
         end;
    until lehet<>1;
   until lehet<>2;
  hiba:;
  if errcode <> 0 then
    writeln(' Error ',errcode,', too many input elements ');
  dispose(xw);
  dispose(yw);
  dispose(xr);
  dispose(yr);
  dispose(x);
  dispose(y);
  dispose(r);
  dispose(t1);
  dispose(t2);
  dispose(ind);
  dispose(indt);
  dispose(time);
  dispose(self);
  dispose(bcells);
  dispose(nonself);
  dispose(marrow);
  dispose(xt);
  dispose(yt);
  dispose(rt);
  dispose(tht1);
  dispose(tht2);
  dispose(abody);
  dispose(thelper);
  dispose(ilns);
end.  {foprogram vege}

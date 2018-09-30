unit G_F;

interface

uses Math,SysUtils;

 type
 bvec= array of byte;
  GF=class
    private
     q: array of integer;
      m: byte;
      value:integer;
      alphapow,logalpha:array of array of integer;

    public
      Constructor Create;   overload;
      constructor Create(qvalue: integer);    overload;
      constructor Create(qvalue, inexp: integer);  overload;
      constructor Create(const ingf: GF ); overload;
      procedure  set_size(qvalue: integer);
      procedure  set_field_value_by_pol(qvalue,inpol : integer);  overload;
      procedure  set_field_value(qvalue, inexp : integer); overload;
      procedure  set_field_value(qvalue : integer;const vectorspace:bvec); overload;
      function  get_vectorspace:bvec;
      procedure  Add(const ingf: GF);
      procedure  Sub(const ingf: GF);
      procedure  Mul(const ingf: GF);
      procedure  Divv(const ingf: GF);

      procedure Copy(const ingf: GF); overload;
      procedure Copy(inexp:integer); overload;
      function  get_value   :integer;
      function  get_size:integer;
      function  ValueToString:string;

  end;

  function Equal(const op1,op2:GF): boolean;
  function NOTEqual(const op1,op2:GF): boolean;
  function Add(const op1,op2:GF):GF;
  function Mul(const op1,op2:GF):GF;
  function Sub(const op1,op2:GF):GF;
  function Divv(const op1,op2:GF):GF;


implementation

constructor GF.Create;
var i:integer;

begin
  SetLength(q,15);

  for i:=0  to 14 do
    q[i]:=1 shl i;

    m:=0;
end;

constructor GF.Create(qvalue: integer);
var i:integer;
begin
  SetLength(q,15);
  m:=0;

  for i:=0  to 14 do
    q[i]:=1 shl i;

  if (qvalue=0) then// qvalue==0 gives the zeroth element
			value:=-1
  else set_size(qvalue);


end;

constructor GF.Create(qvalue, inexp: integer);
var i:integer;
begin
  SetLength(q,15);
  //m:=0;

  for i:=0  to 14 do
    q[i]:=1 shl i;

  set_field_value(qvalue,inexp);


end;

// Copy constructor
constructor GF.Create(const ingf: GF );
var i:integer;

begin
  m:=0;
  SetLength(q,15);

  for i:=0  to 14 do
    q[i]:=1 shl i;

  m:=ingf.m;
  value:=ingf.value;

  set_size(1 shl m);
end;

procedure GF.set_field_value(qvalue,inexp : integer);
begin
		set_size(qvalue);
		value:=inexp;
end;

procedure GF.set_field_value_by_pol(qvalue,inpol : integer);
begin
		set_size(qvalue);
		value:= logalpha[m,inpol];
end;

procedure GF.set_field_value(qvalue: integer;const vectorspace:bvec);
var i,temp,sizebvec:integer;
begin
  temp:=0;
  set_size(qvalue);
  sizebvec:=Length(vectorspace);

	if  (sizebvec<=m) then
  begin
	for i:=0 to sizebvec-1 do
    if (vectorspace[i]=1) then
      temp:=temp +  ( 1 shl (sizebvec-i-1) )   ;

    value:= logalpha[m,temp];
   end;
end;

 function GF.get_vectorspace:bvec;
 var x,i:integer;
     temp:bvec;
 begin
  SetLength(temp,m);

   if (value = -1) then	x:= 0
  else           	x:= alphapow[m,value];

  for i:=m-1 downto 0 do
  begin
	temp[i]:= (x and 1);
	x:= (x shr 1);
  end;

   result:=temp;
 end;

function  GF.get_value :integer;
begin
  result:=value;
end;

procedure GF.set_size(qvalue: integer) ;
var i,mtemp,reduce,n,temp: integer;
  const reducetable: array [0..12] of integer = (3,3,3,5,3,9,29,17,9,5,83,27,43);
begin
  mtemp:=trunc(log10(qvalue)/log10(2.0));
//* Construct gf(q), q=2^m. From Wicker, "Error Control Systems  for digital communication and storage" pp. 463-465 */
	m:=mtemp;
  SetLength(alphapow,m+1);
    SetLength(logalpha,m+1);
  for i := 0 to m do
  begin
   SetLength(alphapow[i],qvalue);
   SetLength(logalpha[i],qvalue);
  end;

  if (m = 1) then // gf(2), special case
  begin
	 alphapow[1,0]:=1;
	 logalpha[1,0]:=-1;
   logalpha[1,1]:=0;
  end
  else
  begin
    reduce:=reducetable[m-2];
	  alphapow[m,0]:=1; // alpha^0 = 1
     logalpha[m,0]:=-1;
	  for n:=1 to (1 shl m)-2 do
    begin
      temp:=alphapow[m,n-1];
	 	  temp:=(temp shl 1); // multiply by alpha
	 	  if ((temp and (1 shl m))<>0) then // contains alpha**m term
      begin
	 	    alphapow[m,n]:=(temp and not(1 shl m)) xor reduce;
      end
      else 
	 	    alphapow[m,n]:=temp; // if no alpha**m term, store as is
		// create table to go in opposite direction
	 	    logalpha[m,n]:=-1; // special case, actually log(0)=-inf

    end;

	  for n:=0 to ((1 shl m)-2) do
      logalpha[m,alphapow[m,n]]:=n;

	end;
 end;

   function GF.get_size:integer;
   begin
   if (m<>0) then  Result:=q[m] else Result:=0;
    end;


procedure GF.Add(const ingf: GF);
begin

if (value = -1) then
begin
	value:=ingf.value;
	m:=ingf.m;
end
else if ((ingf.value <> -1) and (m=ingf.m)) then
  value:=logalpha[m,alphapow[m,value] xor alphapow[m,ingf.value]];

end;

 procedure GF.Mul(const ingf: GF);
 begin
  if ((value = -1) or (ingf.value = -1)) then
	  value:=-1
  else
 if   (ingf.m = m) then
  value:=(value+ingf.value)mod(q[m]-1);
 end;

  procedure GF.Sub(const ingf: GF);
  begin
     Add(ingf);
  end;

  procedure GF.Divv(const ingf: GF);
  begin

  if (ingf.value<>-1) then  

  if (value = -1) then
	value:=-1
  else

  if (ingf.m = m) then
	value:=(value-ingf.value+q[m]-1)mod(q[m]-1);

  end;

procedure GF.Copy(const ingf: GF);
begin
  m:=ingf.m;
  value:=ingf.value;
  set_size(1 shl m);
end;

    procedure GF.Copy(inexp:integer);
    begin

    if ( (m>0) and (inexp>=-1) and ( inexp<(q[m]-1)) ) then
      value:=inexp;

    end;

 function GF.ValueToString:string;
 begin
  Result:='alpha^'+IntToStr(value);
 end;

 function Equal(const op1,op2:GF): boolean;
 begin
 if ((op1.value=-1) and (op2.value=-1)) then
 Result:=true
 else
 if  ((op1.value=op2.value) and (op1.m=op2.m))then
  Result:=true
 else Result:=false;

 end;

  function NOTEqual(const op1,op2:GF): boolean;
  begin
    Result:=not(Equal(op1,op2));
  end;

   function Add(const op1,op2:GF):GF;
   begin
    Result:=GF.Create(op1);
    Result.Add(op2);

   end;

function Mul(const op1,op2:GF):GF;
begin
  Result:=GF.Create(op1);
  Result.Mul(op2);
end;

function Sub(const op1,op2:GF):GF;
begin
 Result:=GF.Create(op1);
 Result.Add(op2);
end;

function Divv(const op1,op2:GF):GF;
begin
  Result:=GF.Create(op1);
  Result.Divv(op2);
end;




end.

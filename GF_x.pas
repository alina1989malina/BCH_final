unit GF_x;

interface
uses G_F, Math,SysUtils;

 type
 bvec= array of byte;
 ivec= array of integer;
  GFX=class
    private
      q, degree : integer;
    public
      coeffs: array of GF;
      Constructor Create;   overload;
      constructor Create(qvalue: integer);  overload;
      constructor Create( qvalue, indegree: integer);    overload;
      constructor Create(qvalue : integer; const invalues : bvec);  overload;
      constructor Create(qvalue : integer; const invalues : string);  overload;
      constructor Create(const ingfx: GFX ); overload;

      procedure set_field_values(qvalue: integer; const invalues: ivec );  overload;
      procedure set_field_values(qvalue: integer; const invalues: string ); overload;
      function get_size :integer;
      function  get_degree:integer ;  // Return degree of gf(q)[x]
      procedure set_degree(indegree:integer);
      function get_true_degree:integer;
      procedure clear;
      procedure Copy(const ingfx: GFX); //overload;
      procedure Add(const ingfx: GFX);
      procedure Sub(const ingfx: GFX);
      procedure Mul(const ingfx: GFX); overload;

      function Count(const ingf:GF):GF;
      function toString : string ;

  end;

      function Add(const op1, op2:GFX):GFX;
      function Sub(const op1, op2:GFX):GFX;
      function Mul(const op1, op2:GFX):GFX;overload;
      function Divv(const c,g: GFX):GFX;overload;
      function Modd(const c,g: GFX):GFX;
      function Mul(const ingfx: GFX; const ingf:GF):GFX;  overload;
      function Divv(const ingfx: GFX; const ingf:GF):GFX;    overload;
implementation

constructor GFX.Create;
begin
  degree:=-1;
  q:=0;
end;

constructor GFX.Create(qvalue: integer);
 begin
 q:=qvalue;
end;

 constructor GFX.Create( qvalue, indegree: integer);
 var i :integer;
 begin
  q:=qvalue;
  SetLength(coeffs,indegree+1);
  degree:=indegree;
  for i:=0 to degree do
	  coeffs[i]:=GF.Create(q,-1);
 end;

 //Copy
constructor GFX.Create(const ingfx: GFX );
  var i,b:integer;
begin

    degree:=ingfx.degree;
    q:=ingfx.q;
    SetLength(coeffs,degree+1);

    for i:=0 to degree do
	   coeffs[i]:=GF.Create(q,ingfx.coeffs[i].get_value);

end;

constructor GFX.Create(qvalue:integer; const invalues: bvec);
var input:ivec;
    i,len  :integer;
begin
  len:=  Length(invalues) ;
  SetLength(input,len);
  for i:=0 to len-1 do
    input[i]:=  invalues[i]   ;
  set_field_values(qvalue,input);
end;

constructor GFX.Create( qvalue:integer; const invalues: string);
begin
  set_field_values(qvalue,invalues);
end;

procedure GFX.set_field_values(qvalue: integer; const invalues: ivec )  ;
var b,i:integer;
begin
  degree:=  Length(invalues)-1;
  Setlength(coeffs,degree+1);
  for i:=0 to degree do
 	  coeffs[i]:=GF.Create(qvalue,invalues[i]);

  q:=qvalue;
end;


procedure GFX.set_field_values(qvalue: integer; const invalues: string ) ;
var  size:byte;
     vec:ivec;
     i:integer;
     forconvert:GF;
begin
  size:=Length(invalues);
  Setlength(vec,size);
  forconvert:=GF.Create(qvalue);

  for i:=1 To size Do
  begin
    forconvert.set_field_value_by_pol(qvalue,StrToInt(invalues[i]));
    vec[i-1]:= forconvert.get_value;
  end;

  set_field_values(qvalue,vec);
end;

   function GFX.get_size :integer;
   begin
     Result:=q;
   end;

     function  GFX.get_degree:integer ;
     begin
     Result:=degree;
   end;

    procedure GFX.set_degree(indegree:integer);
    begin
     SetLength(coeffs,indegree+1);
  degree:=indegree;
  end;

   function GFX.get_true_degree:integer;
   var i:integer;
   begin
   i:=degree;

  while(coeffs[i].get_value=-1) do
  begin
  i:=i-1;
	if (i=-1) then
	  break;
   end;

   Result:=i;

   end;

     procedure GFX.clear;
     var i:integer;
     begin
     for i:=0 to degree do
	   coeffs[i].set_field_value(q,-1);
     end;



  procedure GFX.Copy(const ingfx: GFX);
   var i,Olddegree:integer;
 begin
    Olddegree:=degree;
    degree:=ingfx.degree;
    q:=ingfx.q;
    SetLength(coeffs,degree+1);

    for i:=0 to degree do
      if(i<=Olddegree) then
	   coeffs[i].set_field_value(q,ingfx.coeffs[i].get_value)
    else
	   coeffs[i]:=GF.Create(q,ingfx.coeffs[i].get_value);



    (*for i:=0 to Olddegree-1 do
	   coeffs[i].set_field_value(q,ingfx.coeffs[i].get_value);
    for i:=Olddegree to degree do
	   coeffs[i]:=GF.Create(q,ingfx.coeffs[i].get_value);
      *)
  end;



  procedure GFX.Add(const ingfx: GFX);
  var i,j:integer;
  begin
      if (q = ingfx.q) then
      begin
    if (ingfx.degree > degree) then
    begin

	SetLength(coeffs,ingfx.degree+1);
	// set new coefficients to the zeroth element

	for j:=degree+1 to ingfx.degree do
   coeffs[j]:=GF.create(q,-1);

	degree:=ingfx.degree;
      end;

    for i:=0 to ingfx.degree do
      coeffs[i].Add(ingfx.coeffs[i]);

      end;
  end;

    procedure GFX.Sub(const ingfx: GFX);
    begin
       if (q = ingfx.q) then
     Add(ingfx);
    end;

    procedure GFX.Mul(const ingfx: GFX);
    var i,j:integer;
         tempcoeffs:array of GF;
  begin
     if (q = ingfx.q) then
      begin
      Setlength(tempcoeffs,degree+1);
  for i:=0 to degree do
  tempcoeffs[i]:=GF.Create(coeffs[i]);

    SetLength(coeffs,degree+ingfx.degree+1);
    for j:=0 to degree do
	  coeffs[j].set_field_value(q,-1);     // set coefficients to the zeroth element (log(0)=-Inf=-1)

     for j:=degree+1 to (degree+ingfx.degree) do     // and create new coeffs
	  coeffs[j]:=GF.Create(q,-1);

    for i:=0 to degree do
	for j:=0 to ingfx.degree do
   coeffs[i+j].Add(G_F.Mul(ingfx.coeffs[j],tempcoeffs[i]));


    degree:=degree+ingfx.degree;
      end;
  end;

    function Divv(const c,g: GFX):GFX;
    var q_c,q_g,q,tempdegree,gdegree,degreedif,i:integer;
        temp,m,divisor:GFX;
        Zeroo:GF;
  begin

    q_c:= c.get_size;
     q_g:=g.get_size;
     q:=c.get_size;

     if (q_c=q_g) then
     begin
      Zeroo:=GF.Create(q,-1);
      temp:=GFX.Create(c);
      tempdegree:= temp.get_true_degree;
      gdegree:= g.get_true_degree;

      if (gdegree=-1) then
          Result:=GFX.Create(q_c,0)
      else
      if (tempdegree>=gdegree) then
      begin
        degreedif:= tempdegree - gdegree;
        m:=GFX.Create(q,degreedif);
     divisor:=GFX.Create(q,degreedif);

    for i:=0 to  c.get_degree-1 do
    begin
    	m.coeffs[degreedif]:= G_F.Divv(temp.coeffs[tempdegree],g.coeffs[gdegree]);
    	divisor.set_degree(degreedif);
	    divisor.clear;
	divisor.coeffs[degreedif]:= m.coeffs[degreedif];

 	temp.Sub(Mul(divisor,g));
	tempdegree:= temp.get_true_degree;
	degreedif:= tempdegree - gdegree;
	if ( (degreedif<0) or ((temp.get_true_degree=0) and (Equal(temp.coeffs[0],Zeroo) ) )) then
	    break;

    end;
    Result:=GFX.Create(m);
     end
     else
      Result:=GFX.Create(q_c,0);
    end
     else Result:=GFX.Create;
  end;

    function Modd(const c,g: GFX):GFX;
    var q_c,q_g,q,tempdegree,gdegree,degreedif,i:integer;
        temp,m,divisor:GFX;
        Zeroo:GF;
  begin

     q_c:= c.get_size;
     q_g:=g.get_size;
           q:=q_c;
     if (q_c=q_g) then
     begin

      Zeroo:=GF.Create(q,-1);

      gdegree:= g.get_true_degree;

      if (gdegree=-1) then
          temp:=GFX.Create(q_c,0)
      else
      begin
        temp:=GFX.Create(c);
        tempdegree:= temp.get_true_degree;

        if (tempdegree>=gdegree) then
     begin

     degreedif:= tempdegree - gdegree;
     m:=GFX.Create(q,degreedif);
     divisor:=GFX.Create(q,degreedif);

    for i:=0 to  c.get_degree-1 do
    begin
    	m.coeffs[degreedif]:= G_F.Divv(temp.coeffs[tempdegree],g.coeffs[gdegree]);
    	divisor.set_degree(degreedif);
	    divisor.clear;
	divisor.coeffs[degreedif]:= m.coeffs[degreedif];

 	temp.Sub(Mul(divisor,g));
	tempdegree:= temp.get_true_degree;
	degreedif:= tempdegree - gdegree;
	if ( (degreedif<0) or ((temp.get_true_degree=0) and (Equal(temp.coeffs[0],Zeroo) ) )) then
	    break;

    end;
     end;
      end;
    Result:=GFX.Create(temp);

     end
     else Result:=GFX.Create;  

  end;

  function Mul(const op1,op2:GFX):GFX;
  begin
    Result:=GFX.Create(op1);
          Result.Mul(op2);
  end;



  function GFX.Count(const ingf:GF):GF;
  var i:integer;
      temp,ingfpower:GF;
  begin
  if (q = ingf.get_size) then
  begin
     temp:=GF.Create(coeffs[0]);
      ingfpower:=GF.Create(ingf);

    for i:=1 to degree do
    begin
	temp.Add(G_F.Mul(coeffs[i],ingfpower));
	ingfpower.Mul(ingf);
   end;

     result:= GF.Create(temp);
  end;
     //result:=GF.Create;
  end;


    function GFX.toString : string ;
    var i,terms:integer;
        Zero0,One:GF;
         s:string;
        begin
        terms:=0;

        Zero0:=GF.Create(q,-1);
        One:=GF.Create(q,0);

      for i:=0 to degree do
  	if ( NOTEqual(coeffs[i],Zero0) ) then
    begin
	    if (terms<>0) then Result:=Result+'+' else terms:=terms+1;

	    if (Equal(coeffs[i],One)) then Result:=Result+'x^'+IntToStr(i)
	 	 else Result:=Result+coeffs[i].ValueToString+'*x^'+IntToStr(i);
     end;

    if (terms = 0) then Result:='zero polynom';

        end;


  
  function Add(const op1, op2:GFX):GFX;
  begin
    Result:=GFX.Create(op1);
          Result.Add(op2);
  end;


  
  function Sub(const op1, op2:GFX):GFX;
  begin
      Result:=GFX.Create(op1);
          Result.Add(op2);
  end;


    function Mul(const ingfx:GFX; const  ingf:GF):GFX;
     var temp:GFX;
         i:integer;
     begin
       if (ingf.get_size = ingfx.q) then
       begin
     temp:=GFX.Create(ingfx);

    for i:=0 to ingfx.get_degree do
   	temp.coeffs[i].Mul(ingf);

    Result:=GFX.Create(temp);
       end
       else  Result:=GFX.Create;


     end;


     function Divv(const ingfx:GFX;const ingf:GF):GFX;
     var temp:GFX;
         i:integer;
             begin
       if (ingf.get_size = ingfx.q) then
       begin
      temp:=GFX.Create(ingfx);

    for i:=0 to ingfx.get_degree do
    temp.coeffs[i].Divv(ingf);

    Result:=GFX.Create(temp);
       end
     else  Result:=GFX.Create;
             end;








end.

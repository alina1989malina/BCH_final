unit B_C_H;

interface

uses Math, G_F, GF_x,  SysUtils;

 type
 bvec= array of byte;
  BCH=class
    private
       n,k,t:integer;
       g: GFX;
    public
      Constructor Create (in_n, in_k, in_t: integer;const genpolynom : ivec );  overload  ;
     // constructor Create(in_n, in_k, in_t: integer;var genpolynom : string ); overload  ;
     // Constructor Create (var genpolynom : ivec );  overload  ;
     // Constructor Create (var zeroes:GFX);  overload  ;
      constructor Create(const inbch: BCH);   overload  ;
      function   Showg: string;
      function Encode(const  inuncodedbits : bvec) : bvec;
      function Decode(const  incodedbits : bvec) : string;//bvec;
      function Sindrom(const  validatepol : GFX) : GFX;
     function GetCheckPolinom : GFX;
  end;


implementation

constructor BCH.Create(in_n, in_k, in_t: integer;const  genpolynom : ivec );
var  exponents:ivec;
          temp: bvec;
          i,j:integer;
begin
  n := in_n;
  k := in_k;
  t := in_t;
    //fix the generator polynomial g(x).
  Setlength(exponents,3*Length(genpolynom));

  SetLength(temp,3*Length(genpolynom));
  for i := 0 to Length(genpolynom)-1 do
  begin
    t:=genpolynom[i];
    for j := 2 downto 0 do
      temp[i*3+2-j]:= (t and (1 shl j))shr(j)  ;       //oct2bin
  end;

  for i:=0 to Length(temp)-1 do
 	  exponents[i] := -1 + temp[Length(temp)-i-1];

  g:=GFX.Create(n+1,Length(temp));
  g.set_field_values(n+1,exponents);

end;

function   BCH.Showg: string;
Begin
  Result:=g.tostring;
End;


  (*constructor BCH.Create(in_n, in_k, in_t: integer;var genpolynom : string );
  begin
    n := in_n;
    k := in_k;
    t := in_t;
  //fix the generator polynomil g(x).
    g:=GFX.Create(n+1,genpolynom);
  
  end;  *)


  // Copy constructor
  constructor BCH.Create(const inbch: BCH);
  begin
    n := inbch.n;
    k := inbch.k;
    t := inbch.t;
    g:=GFX.Create(inbch.g);
  end;

  function BCH.Encode(const  inuncodedbits : bvec) : bvec;
  var i,j,itterations,degree, difference :integer;
          m,c: GFX;
          outbit,mbit,cbit, uncodedbits: bvec;
          one:GF;
  begin

  setlength(uncodedbits,length(inuncodedbits));
    for i:=0 to Length(inuncodedbits) - 1 do
      uncodedbits[i]:= inuncodedbits[i];
      
    itterations := ( length(uncodedbits) div k )+1;
    setlength(uncodedbits,itterations*k)  ;
    for i:=length(uncodedbits) to itterations*k-1 do
       uncodedbits[i]:=0;

    m:= GFX.Create(n+1,k);
    c:= GFX.Create(n+1,n);

    SetLength(Result,itterations*n);
    SetLength(mbit,k);
    SetLength(cbit,n);

    One:=GF.Create(n+1,0);

    for i:=0 to itterations-1 do
    begin

  	//Fix the message polynom m(x). mbit = uncodedbits.mid(i*k,k);
      for j:=0  to k - 1 do
        mbit[j]:=uncodedbits[i*k+j];

       for j:=0 to k-1  do
       begin
  	    degree:= -1 + mbit[j];
  	    m.coeffs[j].set_field_value(n+1,degree);
  	    end;

      //Fix the outputbits cbit.
  	     c:=GF_x.Mul(g,m);

  	  for j:=0 to n-1 do
  	    if ( Equal(c.coeffs[j],One) ) then 	cbit[j]:= 1
  	    else cbit[j]:= 0;

  	 for j:=0 to n-1 do
     Result[i*n+j]:=cbit[j];

    end;


  end;


  function BCH.Decode(const  incodedbits : bvec) : string;//bvec ;
  //function BCH.Decode(const  incodedbits : bvec) : string;
  var j, i, degree, kk, foundzeros, cisvalid,itterations:integer;
    r,c,m,S,Lambda,OldLambda, Tpoly, Ohmega, One, x, x2:GFX;
    delta, temp, Zero:GF;
    codedbits, outvec, rbin, rmin ,mbin: bvec     ;
    errorpos : ivec;
  begin

  setlength(codedbits,length(incodedbits));
    for i:=0 to Length(incodedbits) - 1 do
      codedbits[i]:= incodedbits[i];

    itterations := ( length(codedbits) div n )+1;
    setlength(codedbits,itterations*n)  ;
    for i:=length(codedbits) to itterations*n-1 do
       codedbits[i]:=0;

    SetLength(outvec, itterations*k);
    SetLength(rbin, n);
    SetLength(mbin, k);

    r:=GFX.Create(n+1,n-1);
    c:=GFX.Create(n+1,n-1);
    m:=GFX.Create(n+1,k-1);
    S:=GFX.Create(n+1,2*t) ;

    OldLambda:=GFX.Create(n+1,t+1);
    Tpoly:=GFX.Create(n+1);
    Ohmega:=GFX.Create(n+1);
    Zero:=GF.Create(n+1,-1) ;
    One :=GFX.Create(n+1,'1');
    x :=GFX.Create(n+1,'01');
    x2 :=GFX.Create(n+1,'001');

    delta:=GF.Create(n+1);


    for i:=0 to  itterations-1 do
    begin   //Fix the received polynomial r(x)


	    //  rbin := codedbits.mid(i*n,n);
       for j:=0  to n - 1 do
        rbin[j]:=codedbits[i*n+j];

	    for j:=0 to n-1 do
      begin
	    degree := -1 + rbin[j];
	    r.coeffs[j] := GF.Create(n+1,degree);
      end;



	//Fix the syndrome polynomial S(x).
	S.coeffs[0] := GF.Create(n+1,-1);
	for j:=1 to 2*t do
	    S.coeffs[j] :=  r.count(GF.Create(n+1,j));


    (*
	if (S.get_true_degree() >= 1) { //Errors in the received word
	    //Itterate to find Lambda(x).
	    kk = 0;
	    Lambda = gfx(n+1,(char*)//"0");
	    //Tpoly = gfx(n+1,(char*)"0");
	    //while (kk<t) {
		(*Ohmega = Lambda * (S + One);
		delta = Ohmega[2*kk+1];
		OldLambda = Lambda;
		Lambda = OldLambda + delta*( gfx(n+1,(char*)//"-1 0")*Tpoly );
		(*if ((delta == gf(n+1,-1)) || (OldLambda.get_degree() > kk)) {
		//    Tpoly = gfx(n+1,(char*)//"-1 -1 0") * T;
		(*} else {
		 //   Tpoly = ( gfx(n+1,(char*)//"-1 0") * OldLambda ) / delta;
		//}
		//kk = kk + 1;
	  //  }        *)

    Lambda := GFX.Create(n+1,'1');


  if (S.get_true_degree >= 1) then //Errors in the received word
  begin            //Itterate to find Lambda(x).
    kk := 0;
    //Lambda := GFX.Create(n+1,'1');
    Tpoly := GFX.Create(n+1,'1');
    while (kk<t) do
    begin

		  Ohmega := GF_x.Mul(Lambda , GF_x.Add(S,One));
      delta.Copy(Ohmega.coeffs[2*kk+1]);
		  OldLambda.Copy(Lambda);
		  Lambda := GF_x.Add(OldLambda , GF_x.Mul( GF_x.Mul(x,Tpoly), delta ));
		    if ((Equal(delta ,zero)) or (OldLambda.get_degree > kk)) then
		      Tpoly.Mul(x2)
        else
		      Tpoly := GF_x.Divv((GF_x.Mul ( x , OldLambda ) ), delta);

		  kk := kk + 1;
	  end;

      SetLength(errorpos,Lambda.get_true_degree);
        foundzeros := 0;
         for j:=0 to n-1 do
         begin
           temp:=Lambda.count(GF.Create(n+1,j));
           if( Equal(temp,GF.Create(n+1,-1))) then
           begin
             errorpos[foundzeros] :=((n-j) mod n);
             inc(foundzeros) ;
             if (foundzeros >= Lambda.get_true_degree) then break;
           end;
         end;

      for j := 0 to foundzeros - 1 do
        rbin[errorpos[j]]:=(rbin[errorpos[j]]+1) mod 2;


   //end;

	    //Reconstruct the corrected codeword.
	   for j:=0 to n-1 do
     begin
		    degree := -1 + rbin[j];
        c.coeffs[j] := GF.create(n+1,degree);
      end;





	    //Code word validation.

  S.coeffs[0]:=GF.Create(n+1,-1);
  for j := 1 to 2*t do
    S.coeffs[j]:= c.Count(GF.Create(n+1,j));

  if( S.get_true_degree<=0) then   //c(x) is a valid codeword.
    cisvalid:=true
  else  cisvalid:=false;

 end

 else
 begin
   c:=GFX.Create(r);
   cisvalid:=true;
 end;

	//Construct the message bit vector.
	//if (cisvalid) { //c(x) is a valid codeword.
	  //  if (c.get_true_degree() > 1) {
	 //	m = divgfx(c,g);
		//mbin.clear();
		//for (j=0; j<=m.get_true_degree(); j++) {
		//    if ( m[j] == gf(n+1,0) ) {
	 //		mbin(j) = 1;
		//    }
	 //	}
	  //  } else { //The zero word was transmitted
		//mbin = zeros_b(k);
	 //	m = gfx(n+1,(char*)"-1");
	 //   }
	//} else { //Decoder failure.
	 //   mbin = zeros_b(k);
	 //   m = gfx(n+1,(char*)"-1");
	//}
 //	out.replace_mid(i*k,mbin);
  //  }
   // return out;
//}


    if(cisvalid=true) then    //c(x) is a valid codeword.
      if  (c.get_true_degree > 1  ) then
      begin
        m: = divv(c,g);
        mbin.clear;
        for j:=0 to m.get_true_degree do
          if(Equal(m.coeffs[j],GF.Create(n+1,0)))
            mbin.coeffs[j]:=1;
      end
      else  //The zero word was transmitted

Result:=Result+c.toString;
    end;

 end;

 function BCH.GetCheckPolinom : GFX;
 var mn:GFX;
 begin
 mn:=GFX.Create(n+1,n);
 mn.coeffs[0].set_field_value(n+1,0);
  mn.coeffs[n].set_field_value(n+1,0);

  Result:=GFX.Create(n+1);
  Result:=Divv(mn,g);
 end;

  function BCH.Sindrom(const  validatepol : GFX) : GFX;
  begin
     Result:=GFX.Create(n+1);
  Result:=Modd(validatepol,g);
  end;

end.

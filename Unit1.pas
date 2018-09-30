unit Unit1;

interface

uses
  Windows, Messages, SysUtils, Variants, Classes, Graphics, Controls, Forms,
  Dialogs, StdCtrls, G_F, GF_x, B_C_H, Math, ExtCtrls;

type
  TForm1 = class(TForm)
    Label34: TLabel;
    Label35: TLabel;
    Edit7: TEdit;
    Button5: TButton;
    Label36: TLabel;
    Label37: TLabel;
    Label38: TLabel;
    Label39: TLabel;
    Label40: TLabel;
    Label41: TLabel;
    Edit11: TEdit;
    Label42: TLabel;
    Label43: TLabel;
    Edit12: TEdit;
    Button6: TButton;
    Button7: TButton;
    Label44: TLabel;
    Label45: TLabel;
    Label50: TLabel;
    Label51: TLabel;
    Label52: TLabel;
    Label54: TLabel;
    Label55: TLabel;
    Label56: TLabel;
    Label57: TLabel;
    Label1: TLabel;
    Label2: TLabel;
    Label3: TLabel;
    Label4: TLabel;
    Label5: TLabel;
    Label6: TLabel;
    Label7: TLabel;


    procedure Button5Click(Sender: TObject);
    procedure Button6Click(Sender: TObject);
    procedure Button8Click(Sender: TObject);
    procedure Button7Click(Sender: TObject);


  private
    { Private declarations }
  public

    { Public declarations }
  end;

var
  Form1: TForm1;
   GF1,GF2,GF3,GF4,GF5,GF6,GF7,GF8,GF9:GF;
   GFX1,GFX2,GFX3,GFX4,GFX5,GFX6,GFX7,GFX8,GFX9:GFX;
   BCH1, BCH2: BCH;
    qvalue:integer;
implementation

{$R *.dfm}










procedure TForm1.Button5Click(Sender: TObject);
var s,temp, invtemp:string;
size,i,k,n,t:integer;
invec:ivec;
begin
  temp:=Edit7.Text;
  size:=Length( temp);

 invtemp:='';
 for i:=1 To size Do
 begin
   s:= temp[Size-i+1];  //inversion
   invtemp:=invtemp+s;
 end;
 n:=31;
 k:=21;
 t:=2;
  //BCH1:=BCH.Create(n,k,t,invec);
  BCH1:=BCH.Create(invtemp);
 // BCH2:=BCH.Create(BCH1);
 Label55.Caption:=inttostr(BCH1.n);
 Label56.Caption:=inttostr(BCH1.k);
 Label57.Caption:=inttostr(BCH1.t);
   Label4.Caption:=inttostr(BCH1.chain_begin);
     Label5.Caption:=inttostr(2*(BCH1.t));
  Label44.Caption:=BCH1.Showg;
   Label51.Caption:=BCH1.GetCheckPolinom.tostring;

end;

procedure TForm1.Button6Click(Sender: TObject);
var infpolstr,s,codepol, invcodepol:string ;
invec,codevec: bvec;
i,size,k:integer;
infpol:GFX;
begin
  infpolstr:=Edit11.Text;

  size:=Length( infpolstr);
  SetLength(invec,size);
  invcodepol:='';

 for i:=1 To size Do
 begin
   s:= infpolstr[Size-i+1];
   k:= StrToInt(s);
   invec[i-1]:=k;
   invcodepol:=invcodepol+s;
 end;

 infpol:=GFX.Create(2,invcodepol);
  Label6.Caption:= infpol.toString;

 codevec:=BCH1.Encode(invec);

 for i:=0 To Length(codevec)-1 Do
 begin
   codepol:=codepol+inttostr(codevec[Length(codevec)-i-1]);
 end;

 Label45.Caption:=codepol;
 Edit12.Text:=codepol;

end;

procedure TForm1.Button7Click(Sender: TObject);
var codev,s,s2, invcodev: string  ;
cod, decodevec: bvec;
codevec:GFX;
size,k,i:integer;
codepol:GFX;
begin
  codev:=Edit12.Text;

  size:=Length( codev);
  SetLength(cod,size);
  invcodev:='';

 for i:=1 To size Do
 begin
   s:= codev[Size-i+1];
   k:= StrToInt(s);
   invcodev:=invcodev+s;
   cod[i-1]:=k;
 end;

 codepol:=GFX.Create(2, invcodev);
  Label7.Caption:=codepol.toString;

   //codevec:=GFX.Create(32,codev);
   //Label50.Caption:=BCH1.Sindrom(codevec).toString;

   decodevec:=BCH1.Decode(cod);

   s2:='';
   for i:=0 To Length(decodevec)-1 Do
    s2:= s2+inttostr(decodevec[Length(decodevec)-i-1]);

   Label54.Caption:=s2;


 // :=Edit12.Text;
  //Label50.Caption:=




end;

procedure TForm1.Button8Click(Sender: TObject);
var u0:string;
kk,i,t0,n:integer;
a:array of byte;
 S, Lambda, OldLambda, T, Ohmega,x, x2,one, incS:GFX;
 delta,zero0:GF;
begin
  (*n:=1;
  t0:=3;

  OldLambda:=GFX.Create(n+1,t0+1);
  T:=GFX.Create(n+1);
  Ohmega:=GFX.Create(n+1);
  zero0:=GF.Create(n+1,-1) ;
  One :=GFX.Create(n+1,'1');
  x :=GFX.Create(n+1,'01');
  x2 :=GFX.Create(n+1,'001');

  delta:=GF.Create(n+1);

  u0:=Edit13.Text;
  S:=GFX.Create(n+1,u0);

  incS:= GFX.Create(n+1,Length(u0)+1); *)













	 //   Lambda = gfx(n+1,(char*)"0");
	  //  T = gfx(n+1,(char*)"0");
	  //  while (kk<t) {
	 //	Ohmega = Lambda * (S + One);
	 //	delta = Ohmega[2*kk+1];
	//	OldLambda = Lambda;
	//	Lambda = OldLambda + delta*( gfx(n+1,(char*)"-1 0")*T );
 //		if ((delta == gf(n+1,-1)) || (OldLambda.get_degree() > kk)) {
	//	    T = gfx(n+1,(char*)"-1 -1 0") * T;
 //		} else {
 //		    T = ( gfx(n+1,(char*)"-1 0") * OldLambda ) / delta;
 //		} 



 (*	if (S.get_true_degree >= 1) then //Errors in the received word
  begin            //Itterate to find Lambda(x).
    kk := 0;
    Lambda := GFX.Create(n+1,'1');
    T := GFX.Create(n+1,'1');
    while (kk<t0) do
    begin

		  Ohmega := GF_x.Mul(Lambda , GF_x.Add(S,One));
      delta.Copy(Ohmega.coeffs[2*kk+1]);
		  OldLambda.Copy(Lambda);
		  Lambda := GF_x.Add(OldLambda , GF_x.Mul( GF_x.Mul(x,T), delta ));
		    if ((Equal(delta ,zero0)) or (OldLambda.get_degree > kk)) then
		      T.Mul(x2)
        else
		      T := GF_x.Divv((GF_x.Mul ( x , OldLambda ) ), delta);

		  kk := kk + 1;
	  end;
  end;
  Label49.Caption:=Lambda.toString;    *)


end;

end.

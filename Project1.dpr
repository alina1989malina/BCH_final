program Project1;

uses
  Forms,
  Unit1 in 'Unit1.pas' {Form1},
  G_F in 'G_F.pas',
  GF_x in 'GF_x.pas',
  B_C_H in 'B_C_H.pas';

{$R *.res}

begin
  Application.Initialize;
  Application.MainFormOnTaskbar := True;
  Application.CreateForm(TForm1, Form1);
  Application.Run;
end.

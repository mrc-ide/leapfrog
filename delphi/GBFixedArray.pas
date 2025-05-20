unit GBFixedArray;

interface

uses
  System.SysUtils, System.Generics.Defaults, System.Generics.Collections;

// This is a runtime fixed length array compared to the Delphi built-in
// fixed-length array which must be known at compile time.
// Note this should only be used when you need the data to be stored in
// contiguous memory i.e. for passing by pointer to a DLL.
// Delphi dynamic arrays are jagged. The inner most dimension is contiguous
// but the outer dimensions are arrays of pointers.
type
  TGBFixedArray<T> = Class
  type
    P = ^T;
  strict private
    FData: TArray<T>;
    FDims: TArray<Integer>;
    Size: Integer;
    function GetFlatIndex(Indices: array of Integer): Integer;
    procedure SetValue(Indices: array of Integer; const Value: T);
    function GetValue(Indices: array of Integer): T;
  public
    constructor Create(const Dims: TArray<Integer>);

    property Values[Indices: array of Integer]: T read GetValue write SetValue; default;

    function GetLength(): Integer;

    function DataPointer: P;
    property Data: P read DataPointer;
  end;

implementation

constructor TGBFixedArray<T>.Create(const Dims: TArray<Integer>);
var
  i: Integer;
begin
  FDims := Copy(Dims);
  Size := 1;
  for i := 0 to High(Dims) do
    Size := Size * Dims[i];
  SetLength(FData, Size);
end;

// Column-major storage order as this is what Eigen::Tensor uses by default
function TGBFixedArray<T>.GetFlatIndex(Indices: array of Integer): Integer;
var
  i, Stride, Offset: Integer;
begin
  if Length(Indices) <> Length(FDims) then
    raise Exception.Create('Incorrect number of indices');

  Result := 0;
  Stride := 1;

  for i := 0 to High(FDims) do
  begin
    if (Indices[i] < 0) or (Indices[i] >= FDims[i]) then
      raise Exception.CreateFmt('Index %d out of bounds for dimension %d', [Indices[i], i]);

    Result := Result + (Indices[i] * Stride);
    Stride := Stride * FDims[i];
  end;
end;

procedure TGBFixedArray<T>.SetValue(Indices: array of Integer; const Value: T);
begin
  FData[GetFlatIndex(Indices)] := Value;
end;

function TGBFixedArray<T>.GetValue(Indices: array of Integer): T;
begin
  Result := FData[GetFlatIndex(Indices)];
end;

function TGBFixedArray<T>.GetLength(): Integer;
begin
  Result := Size;
end;

function TGBFixedArray<T>.DataPointer: P;
begin
  if Length(FData) = 0 then
    Result := nil
  else
    Result := @FData[0];
end;

end.

namespace SystemOfEquations
{
    internal class GausMethod
    {
        public uint ColumCount;
        public uint RowCount;

        public GausMethod(uint Row, uint Colum)
        {
            RightPart = new double[Row];
            Answer = new double[Row];
            Matrix = new double[Row][];
            for (var i = 0; i < Row; i++)
                Matrix[i] = new double[Colum];
            RowCount = Row;
            ColumCount = Colum;


            for (var i = 0; i < Row; i++)
            {
                Answer[i] = 0;
                RightPart[i] = 0;
                for (var j = 0; j < Colum; j++)
                    Matrix[i][j] = 0;
            }
        }

        public double[][] Matrix { get; set; }
        public double[] RightPart { get; set; }
        public double[] Answer { get; set; }

        private void SortRows(int SortIndex)
        {
            var MaxElement = Matrix[SortIndex][SortIndex];
            var MaxElementIndex = SortIndex;
            for (var i = SortIndex + 1; i < RowCount; i++)
            {
                if (!(Matrix[i][SortIndex] > MaxElement)) continue;
                MaxElement = Matrix[i][SortIndex];
                MaxElementIndex = i;
            }


            if (MaxElementIndex <= SortIndex) return;

            var Temp = RightPart[MaxElementIndex];
            RightPart[MaxElementIndex] = RightPart[SortIndex];
            RightPart[SortIndex] = Temp;

            for (var i = 0; i < ColumCount; i++)
            {
                Temp = Matrix[MaxElementIndex][i];
                Matrix[MaxElementIndex][i] = Matrix[SortIndex][i];
                Matrix[SortIndex][i] = Temp;
            }
        }

        public int SolveParallel()
        {
            if (RowCount != ColumCount)
                return 1;


            for (var i = 0; i < RowCount - 1; i++)
            {
                SortRows(i);
                for (var j = i + 1; j < RowCount; j++)
                {
                    if (Matrix[i][i] == 0.0) continue;

                    var MultElement = Matrix[j][i]/Matrix[i][i];
                    for (var k = i; k < ColumCount; k++)
                        Matrix[j][k] -= Matrix[i][k]*MultElement;
                    RightPart[j] -= RightPart[i]*MultElement;
                }
            }


            for (var i = (int) (RowCount - 1); i >= 0; i--)
            {
                Answer[i] = RightPart[i];

                for (var j = (int) (RowCount - 1); j > i; j--)
                    Answer[i] -= Matrix[i][j]*Answer[j];

                if (Matrix[i][i] == 0)
                    return RightPart[i] == 0 ? 2 : 1;

                Answer[i] /= Matrix[i][i];
            }
            return 0;
        }

        public int SolveMatrix()
        {
            if (RowCount != ColumCount)
                return 1;

            for (var i = 0; i < RowCount - 1; i++)
            {
                SortRows(i);
                for (var j = i + 1; j < RowCount; j++)
                {
                    if (Matrix[i][i] == 0) continue;
                    var MultElement = Matrix[j][i]/Matrix[i][i];
                    for (var k = i; k < ColumCount; k++)
                        Matrix[j][k] -= Matrix[i][k]*MultElement;
                    RightPart[j] -= RightPart[i]*MultElement;
                }
            }


            for (var i = (int) (RowCount - 1); i >= 0; i--)
            {
                Answer[i] = RightPart[i];

                for (var j = (int) (RowCount - 1); j > i; j--)
                    Answer[i] -= Matrix[i][j]*Answer[j];

                if (Matrix[i][i] == 0)
                    return RightPart[i] == 0 ? 2 : 1;

                Answer[i] /= Matrix[i][i];
            }
            return 0;
        }


        public override string ToString()
        {
            var S = "";
            for (var i = 0; i < RowCount; i++)
            {
                S += "\r\n";
                for (var j = 0; j < ColumCount; j++)
                {
                    S += Matrix[i][j].ToString("F04") + "\t";
                }

                S += "\t" + Answer[i].ToString("F08");
                S += "\t" + RightPart[i].ToString("F04");
            }
            return S;
        }
    }
}
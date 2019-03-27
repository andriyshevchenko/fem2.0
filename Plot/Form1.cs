using FEM;
using MathNet.Numerics.Integration;
using System;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Windows.Forms;
using System.Windows.Forms.DataVisualization.Charting;
using static System.Collections.Generic.Create;
using static System.Collections.Generic.SeriesCreate;
using static System.Linq.Enumerable;
using static System.Math;

namespace Plot
{
    public partial class Form1 : Form
    {
        private double integrate_d_dx_2(int i)
        {
            double a = x[i];
            double b = x[i + 1];
            double h = b - a;
            double a1 = x[i] + x[i + 1];
            double antiderivative(double val) =>
                (4.0*Pow(val,3))/3.0 - val*val*2*a1 +a1*a1*val;
            return (16.0 * (antiderivative(b) - antiderivative(a))) / (h * h);
        }

        private DiffusionConvectionReaction Task =
            new DiffusionConvectionReaction(
                x,
                new DiffusionConvectionReaction.BoundaryCondition(u0: 0, u1: 0),
                f: x => Pow(Cos(PI * x), 2) + 0.005 * PI * PI * Cos(2 * PI * x),
                mu: 0.0025, beta: 0, omega: 1.0, alpha: 1000.0
            );

        public Form1()
        {
            InitializeComponent();
            Task.Solve();
            Task.Calc_Eh();
            chart1.Series.Add(new Series
            {
                ChartType = SeriesChartType.Line,
                Color = Color.Green,
                MarkerStyle = MarkerStyle.Circle,
                MarkerSize = 4,
                BorderWidth = 1,
                Name = "Uh(x)"
            });

            chart1.Series.Add(new Series
            {
                ChartType = SeriesChartType.Line,
                Color = Color.Red,
                MarkerStyle = MarkerStyle.None,
                MarkerSize = 4,
                BorderWidth = 1,
                Name = "Eh(x)",
                ToolTip = "Eh(x)",
            });
            chart1.Series.Add(new Series
            {
                ChartType = SeriesChartType.Line,
                Color = Color.Black,
                BorderWidth = 1,
                Name = "x = 0"
            });
            chart1.Series[3].Points.DataBindXY(array(0.0, 1), array(0.0, 0.0));

            chart1.Series.Add(new Series
            {
                ChartType = SeriesChartType.Line,
                Color = Color.Gray,
                BorderWidth = 1,
                Name = ""
            });

            chart1.ChartAreas[0].AxisX.Minimum = 0;
            chart1.ChartAreas[0].AxisX.Maximum = 1.0;
            chart1.ChartAreas[0].AxisY.Maximum = 1.0;
            chart1.ChartAreas[0].AxisY.Minimum = -0.4;
            chart1.ChartAreas[0].AxisX.MajorGrid.LineWidth = 0;
            chart1.ChartAreas[0].AxisY.MajorGrid.LineWidth = 0;

            Plot();

            Console.WriteLine(GaussLegendreRule.Integrate(k => B_BasisFunc.d_dx(k, 1, x)*B_BasisFunc.d_dx(k, 1, x), 0, 1, 5));
            Console.WriteLine(integrate_d_dx_2(2));
            //
        }

        private void Plot()
        {
            chart1.Series[1].Points.DataBindXY(Task.Elements, array(Task.Elements.Select(x => Task.U(x))));
            chart1.Series[1].ToolTip = "Uh(x), Elements: " + Task.Elements.Count; 
            chart1.Series[2].Points.DataBindXY(error_x_values, error_x_values.Select(x => Task.Error(x)).ToList());
        }

        private static double[] error_x_values = Series(0.0, 1.0, 240);
        private static double[] x = Series(0.0, 1.0, 4);

        private void Chart1_Click(object sender, System.EventArgs e)
        {
            this.Task.StartAdaptationAlgorithm(0.2);
            Plot();
        }
    }
}

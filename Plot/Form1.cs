using System.Data;
using System.Linq;
using System.Windows.Forms;
using static System.Linq.Enumerable;
using static System.Collections.Generic.Create;
using static System.Collections.Generic.SeriesCreate;
using static System.Math;
using FEM;
using System.Windows.Forms.DataVisualization.Charting;
using System.Drawing;

namespace Plot
{
    public partial class Form1 : Form
    {
        DiffusionConvectionReaction Task =
            new DiffusionConvectionReaction(
                _values,
                new DiffusionConvectionReaction.BoundaryCondition(u0: 0, u1: 0),
                f: x => Pow(Cos(PI * x), 2) + 0.005 * PI * PI * Cos(2 * PI * x),
                mu: 0.0025, beta: 0, omega: 1.0, alpha: 1000.0
            );

        public Form1()
        {
            InitializeComponent();

            chart1.Series.Add(new Series
            {
                ChartType = SeriesChartType.Line,
                Color = Color.Green,
                MarkerStyle = MarkerStyle.Circle,
                MarkerSize = 4,
                BorderWidth = 1
            });

            chart1.Series.Add(new Series
            {
                ChartType = SeriesChartType.Line,
                Color = Color.Red,
                MarkerStyle = MarkerStyle.None,
                MarkerSize = 4,
                BorderWidth = 1
            });
            chart1.Series.Add(new Series
            {
                ChartType = SeriesChartType.Line,
                Color = Color.Black,  
                BorderWidth = 1,  
            });

            chart1.ChartAreas[0].AxisX.Minimum = 0;
            chart1.ChartAreas[0].AxisX.Maximum = 1.0;
            chart1.ChartAreas[0].AxisY.Maximum = 1.0;
            chart1.ChartAreas[0].AxisY.Minimum = -0.4;

            chart1.Series[1].Points.DataBindXY(_values, array(_values.Select(x => Task.U(x))));
            double[] xValue = Series(0.0, 1.0, 120);
            Task.Calc_Eh();
            double[] error_xasis = Series(0.0, 1.0, 240);
            chart1.Series[2].Points.DataBindXY(error_xasis, error_xasis.Select(x => Task.ErrorEstimate(x)).ToList());
            chart1.Series[3].Label = "x = 0";
            chart1.Series[3].Points.DataBindXY(array(0.0, 1), array(0.0, 0.0));
        }

        static double[] _values = Series(0.0, 1.0, 6);
    }
}

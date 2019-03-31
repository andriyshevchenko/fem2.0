using FEM;
using MathNet.Numerics.Integration;
using System;
using System.Collections.Generic;
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
        private DiffusionConvectionReaction Task =
            new DiffusionConvectionReaction(
                x,
                new DiffusionConvectionReaction.BoundaryCondition(u0: 0, u1: 0),
                f: x => 100,
                mu: 1, beta: 100, sigma: 0, alpha: 1000.0
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
                MarkerSize = 10,
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
                Color = Color.Gray,
                BorderWidth = 1,
                Name = "Eta(element)",
                ToolTip = "Eta(element)",
                MarkerStyle = MarkerStyle.Circle,
                MarkerSize = 4
            });

            chart1.ChartAreas[0].AxisX.Minimum = 0;
            chart1.ChartAreas[0].AxisX.Maximum = 1.0;
            chart1.ChartAreas[0].AxisY.Maximum = 1.0;
            chart1.ChartAreas[0].AxisY.Minimum = -0.4;
            chart1.ChartAreas[0].AxisX.MajorGrid.LineWidth = 0;
            chart1.ChartAreas[0].AxisY.MajorGrid.LineWidth = 0;

            Plot();

            Console.WriteLine(CourantFunction.integate_fx_fx(1, 0, x));
        }

        private void Plot()
        {
            chart1.Series[1].Points.DataBindXY(Task.Elements, array(Task.Elements.Select(x => Task.U(x))));
            chart1.Series[1].ToolTip = "Uh(x), Elements: " + Task.Elements.Count; 
            //chart1.Series[2].Points.DataBindXY(error_x_values, error_x_values.Select(x => Task.U(x)+Task.Error(x)).ToList());
        }

        private static double[] error_x_values = Series(0.0, 1.0, 240);
        private static double[] x = Series(0.0, 1.0, 5);

        private void Chart1_Click(object sender, System.EventArgs e)
        {
            List<double> eta = this.Task.StartAdaptationAlgorithm(10);
            Plot();
            //chart1.Series[4].Points.DataBindXY(Series(0.0, 1, eta.Count), eta);
        }
    }
}

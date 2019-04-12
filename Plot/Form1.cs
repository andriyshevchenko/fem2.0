using FEM;
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
        private DiffusionConvectionReaction Task2  =
            new DiffusionConvectionReaction(
                x,
                new DiffusionConvectionReaction.BoundaryCondition(u_: 0, u1: 0),
                f: x => Pow(Cos(PI*x),2) + (0.005*PI*PI*Cos(2*PI*x)),
                mu: 0.0025, beta: 0, sigma: 1.0, alpha: 1000.0
            );

        private DiffusionConvectionReaction Task  =
            new DiffusionConvectionReaction(
                x,
                new DiffusionConvectionReaction.BoundaryCondition(u_: 0, u1: 0),
                f: x => 100,
                mu: 1, beta: 100, sigma: 0, alpha: 1000.0
            );

        public Form1()
        {
            report.Controls.Add(table);
            report.Show();
            this.WindowState = FormWindowState.Minimized;

            InitializeComponent();

            Task.Solve();
            Task.Calc_Eh();
            List<string> outputdata = Task.CalcEta();

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
                ChartType = SeriesChartType.Column,
                Color = Color.Gray,
                BorderWidth = 1,
                Name = "Eta(element)",
                ToolTip = "Eta(element)",
                MarkerStyle = MarkerStyle.None,
                MarkerSize = 4
            });
            chart1.Series[4]["PixelPointWidth"] = "4";
            chart1.ChartAreas[0].AxisX.Minimum = 0;
            chart1.ChartAreas[0].AxisX.Maximum = 1.0;
            chart1.ChartAreas[0].AxisY.Maximum = 1.0;
            chart1.ChartAreas[0].AxisY.Minimum = -0.4;
            chart1.ChartAreas[0].AxisX.MajorGrid.LineWidth = 0;
            chart1.ChartAreas[0].AxisY.MajorGrid.LineWidth = 0;

            Plot();

            this.Load += (sender, args) =>
            {
                this.table.Rows.Add(OutputData(outputdata));
            };
        }

        private object[] OutputData(List<string> outputdata)
        {
            Bitmap plot = new Bitmap(Width, Height);
            chart1.DrawToBitmap(plot, new Rectangle(0, 0, Width, Height));
            object[] values = new object[7];
            outputdata.ToArray().CopyTo(values, 0);
            values[values.Length-1] = new Bitmap(plot, new Size((int)(Width * 0.4), (int)(Height * 0.4)));
            return values;
        }

        private void Plot()
        {
            chart1.Series[1].Points.DataBindXY(Task.Elements, array(Task.Elements.Select(x => Task.U(x))));
            chart1.Series[1].ToolTip = "Uh(x), Elements: " + Task.Elements.Count; 
            chart1.Series[2].Points.DataBindXY(error_x_values, error_x_values.Select(x => Task.Error(x)).ToList());
            chart1.Series[4].Points.DataBindXY(Task.Elements.Take(Task.Elements.Count - 1).ToList(), Task.Eta.Select(item => item / 100.0).ToList());
        }

        private static double[] error_x_values = Series(0.0, 1.0, 2000);
        private static double[] x = Series(0.0, 1.0, 3);

        private async void Chart1_Click(object sender, System.EventArgs e)
        {
            await this.Task.StartAdaptationAlgorithm(1, (list) =>
            {
                Plot();
                this.table.Rows.Add(OutputData(list));
                return System.Threading.Tasks.Task.Delay(1000);
            });
            //Plot();
        }
    }
}

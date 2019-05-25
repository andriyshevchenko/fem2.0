using System;
using System.Drawing;
using System.Windows.Forms;

namespace Plot
{
    partial class Form1
    {
        /// <summary>
        /// Required designer variable.
        /// </summary>
        private System.ComponentModel.IContainer components = null;
        private static readonly Color color = Color.FromArgb(246, 247, 250);
        private Form report = new Form()
        {
            WindowState = FormWindowState.Maximized,
            BackColor = color, 
            Padding = new Padding(20,0,0,0)
        };

        private DataGridView table = new DataGridView()
        {
            Font = new Font("Segoe UI", 10, FontStyle.Bold),
            ColumnHeadersBorderStyle = DataGridViewHeaderBorderStyle.None,
            SelectionMode = DataGridViewSelectionMode.CellSelect,
            AllowUserToAddRows = false,
            AllowUserToDeleteRows = false,
            AllowUserToOrderColumns = true, 
            MultiSelect = false,
            AutoSizeRowsMode = DataGridViewAutoSizeRowsMode.AllCells,
            AllowUserToResizeColumns = true,
            ColumnHeadersHeightSizeMode = DataGridViewColumnHeadersHeightSizeMode.DisableResizing,
            AllowUserToResizeRows = false,
            RowHeadersBorderStyle = DataGridViewHeaderBorderStyle.None,
            RowHeadersWidthSizeMode = DataGridViewRowHeadersWidthSizeMode.DisableResizing,
            EditMode = DataGridViewEditMode.EditProgrammatically,
            ReadOnly = true,
            Dock = DockStyle.Fill,
            BorderStyle = BorderStyle.None,
            BackgroundColor = color,
            ForeColor = Color.Black,
            GridColor = color,
            RowHeadersVisible = false,
            CellBorderStyle = DataGridViewCellBorderStyle.Single
        };

        /// <summary>
        /// Clean up any resources being used.
        /// </summary>
        /// <param name="disposing">true if managed resources should be disposed; otherwise, false.</param>
        protected override void Dispose(bool disposing)
        {
            if (disposing && (components != null))
            {
                components.Dispose();
            }
            base.Dispose(disposing);
        }

        #region Windows Form Designer generated code

        /// <summary>
        /// Required method for Designer support - do not modify
        /// the contents of this method with the code editor.
        /// </summary>
        private void InitializeComponent()
        {
            //
            // table
            //
            table.Columns.Clear();
            table.Columns.AddRange(new DataGridViewColumn[]
            {
                new DataGridViewTextBoxColumn()
                {
                    HeaderText = "Elements",
                    AutoSizeMode = DataGridViewAutoSizeColumnMode.AllCells,
                },
                new DataGridViewTextBoxColumn()
                {
                    HeaderText = "New el",
                    AutoSizeMode = DataGridViewAutoSizeColumnMode.AllCells,
                },
                new DataGridViewTextBoxColumn()
                {
                    HeaderText = "||Eh||",
                    AutoSizeMode = DataGridViewAutoSizeColumnMode.AllCells
                },
                new DataGridViewTextBoxColumn()
                {
                    HeaderText = "||Uh||",
                    AutoSizeMode = DataGridViewAutoSizeColumnMode.AllCells
                },
                new DataGridViewTextBoxColumn()
                {
                    HeaderText = "Max η, %",
                    AutoSizeMode = DataGridViewAutoSizeColumnMode.AllCells
                },
                new DataGridViewTextBoxColumn()
                {
                    HeaderText = "P",
                    AutoSizeMode = DataGridViewAutoSizeColumnMode.AllCellsExceptHeader
                },
                new DataGridViewTextBoxColumn()
                {
                    HeaderText = "Integral",
                    AutoSizeMode = DataGridViewAutoSizeColumnMode.AllCells
                },
                new DataGridViewImageColumn()
                {
                    AutoSizeMode = DataGridViewAutoSizeColumnMode.Fill
                }
            });
            table.RowsAdded += (sender, args) =>
            {
                foreach (DataGridViewCell item in table.Rows[args.RowIndex].Cells)
                {
                    item.Style.BackColor = color;
                }

                table.FirstDisplayedScrollingRowIndex = table.RowCount - 1;
            };
            table.SelectionChanged += (object sender, EventArgs args) =>
            {
                table.ClearSelection(); 
            };

            System.Windows.Forms.DataVisualization.Charting.ChartArea chartArea1 = new System.Windows.Forms.DataVisualization.Charting.ChartArea();
            System.Windows.Forms.DataVisualization.Charting.Legend legend1 = new System.Windows.Forms.DataVisualization.Charting.Legend();
            System.Windows.Forms.DataVisualization.Charting.Series series1 = new System.Windows.Forms.DataVisualization.Charting.Series();
            this.chart1 = new System.Windows.Forms.DataVisualization.Charting.Chart();
            ((System.ComponentModel.ISupportInitialize)(this.chart1)).BeginInit();
            this.SuspendLayout();
            // 
            // chart1
            // 
            chartArea1.Name = "ChartArea1";
            this.chart1.ChartAreas.Add(chartArea1);
            this.chart1.Dock = System.Windows.Forms.DockStyle.Fill;
            legend1.Name = "Legend1";
            this.chart1.Legends.Add(legend1);
            this.chart1.Location = new System.Drawing.Point(0, 0);
            this.chart1.Name = "chart1";
            series1.ChartArea = "ChartArea1";
            series1.Legend = "Legend1";
            series1.Name = "Series1";
            this.chart1.Series.Add(series1);
            this.chart1.Size = new System.Drawing.Size(523, 335);
            this.chart1.TabIndex = 1;
            this.chart1.Text = "chart1";
            this.chart1.Click += new System.EventHandler(this.Chart1_Click);
            // 
            // Form1
            // 
            this.AutoScaleDimensions = new System.Drawing.SizeF(6F, 13F);
            this.AutoScaleMode = System.Windows.Forms.AutoScaleMode.Font;
            this.ClientSize = new System.Drawing.Size(523, 335);
            this.Controls.Add(this.chart1);
            this.Name = "Form1";
            this.Text = "Form1";
            this.WindowState = System.Windows.Forms.FormWindowState.Maximized;
            ((System.ComponentModel.ISupportInitialize)(this.chart1)).EndInit();
            this.ResumeLayout(false);

        }

        #endregion

        private System.Windows.Forms.DataVisualization.Charting.Chart chart1;
    }
}
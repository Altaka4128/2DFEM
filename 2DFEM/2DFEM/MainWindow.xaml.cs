using MathNet.Numerics.LinearAlgebra.Double;
using Microsoft.Win32;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Data;
using System.Windows.Documents;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using System.Windows.Navigation;
using System.Windows.Shapes;

namespace _2DFEM
{
    /// <summary>
    /// MainWindow.xaml の相互作用ロジック
    /// </summary>
    public partial class MainWindow : Window
    {
        private FEM fem;

        public MainWindow()
        {
            InitializeComponent();
            Loaded += loadedEvent;
        }

        private void loadedEvent(object sender, RoutedEventArgs e)
        {
            ClearClicked(null, null);
        }

        private void FileOpen1dClicked(object sender, RoutedEventArgs e)
        {
            // 初期化
            ClearClicked(null, null);

            // 入力用のCSVファイルを読み込む
            var diag = new OpenFileDialog();
            diag.Filter = "CSVファイル (*.csv)|*.csv";
            if (diag.ShowDialog() == true)
            {
                FEMData2d femData = new FEMData2d();
                femData.ReadCSVFile1d(diag.FileName);

                fem = new Shell1dFEM(femData);

                // モデルを描画する
                var elems = fem.Elems;
                for (int i = 0; i < elems.Count; i++)
                {
                    var elem = elems[i].ShallowCopy();   // 要素をコピーする

                    // 描画する
                    DrawShell1dElement((Shell1dElement)elems[i], Brushes.LightGreen, Brushes.Blue, 1.0);
                }
            }
        }

        private void FileOpen2dClicked(object sender, RoutedEventArgs e)
        {
            // 初期化
            ClearClicked(null, null);

            // 入力用のCSVファイルを読み込む
            var diag = new OpenFileDialog();
            diag.Filter = "CSVファイル (*.csv)|*.csv";
            if (diag.ShowDialog() == true)
            {
                FEMData2d femData = new FEMData2d();
                femData.ReadCSVFile2d(diag.FileName);

                fem = new Shell2dFEM(femData);

                // モデルを描画する
                var elems = fem.Elems;
                for (int i = 0; i < elems.Count; i++)
                {
                    var elem = elems[i].ShallowCopy();   // 要素をコピーする

                    // 描画する
                    DrawShell2dElement((Shell2dElement)elems[i], Brushes.LightGreen, Brushes.Blue, 1.0);
                }
            }
        }

        private void Analysis1dClicked(object sender, RoutedEventArgs e)
        {
            // 例外処理
            if (fem == null)
            {
                return;
            }

            // FEM解析を実行する
            fem.Analysis();

            // モデルを描画する
            var elems = fem.Elems;
            DenseVector dispElemVector = DenseVector.Create(8, 0.0);
            DenseVector dispVector = fem.DispVector;
            for (int i = 0; i < elems.Count; i++)
            {
                for (int j = 0; j < 4; j++)
                {
                    dispElemVector[2 * j] = dispVector[2 * (elems[i].Nodes[j].No - 1)];
                    dispElemVector[2 * j + 1] = dispVector[2 * (elems[i].Nodes[j].No - 1) + 1];
                }

                // 変形後の要素を計算する
                var elem = elems[i].ShallowCopy();   // 変形前の要素をコピーする
                for (int j = 0; j < 4; j++)
                {
                    elem.Nodes[j].Point.X += dispElemVector[2 * j];
                    elem.Nodes[j].Point.Y += dispElemVector[2 * j + 1];
                }

                // 変形後の要素を描画する
                DrawShell1dElement((Shell1dElement)elem, Brushes.LightCoral, Brushes.Red, 0.5);
            }
        }

        private void Analysis2dClicked(object sender, RoutedEventArgs e)
        {
            // 例外処理
            if (fem == null)
            {
                return;
            }

            // FEM解析を実行する
            fem.Analysis();

            // モデルを描画する
            var elems = fem.Elems;
            DenseVector dispElemVector = DenseVector.Create(16, 0.0);
            DenseVector dispVector = fem.DispVector;
            for (int i = 0; i < elems.Count; i++)
            {
                for (int j = 0; j < 8; j++)
                {
                    dispElemVector[2 * j] = dispVector[2 * (elems[i].Nodes[j].No - 1)];
                    dispElemVector[2 * j + 1] = dispVector[2 * (elems[i].Nodes[j].No - 1) + 1];
                }

                // 変形後の要素を計算する
                var elem = elems[i].ShallowCopy();   // 変形前の要素をコピーする
                for (int j = 0; j < 8; j++)
                {
                    elem.Nodes[j].Point.X += dispElemVector[2 * j];
                    elem.Nodes[j].Point.Y += dispElemVector[2 * j + 1];
                }

                // 変形後の要素を描画する
                DrawShell2dElement((Shell2dElement)elem, Brushes.LightCoral, Brushes.Red, 0.5);
            }
        }

        private void ClearClicked(object sender, RoutedEventArgs e)
        {
            fem = null;
            this.Canvas.Children.Clear();
            DrawXYArrow();
        }

        private void DrawXYArrow()
        {
            // X軸を描画する
            Line xArrow = new Line();
            xArrow.Stroke = Brushes.Red;
            xArrow.StrokeThickness = 4;
            xArrow.StrokeEndLineCap = PenLineCap.Triangle;
            xArrow.X1 = 0;
            xArrow.Y1 = 0;
            xArrow.X2 = 50;
            xArrow.Y2 = 0;
            this.Canvas.Children.Add(xArrow);

            // Y軸を描画する
            Line yArrow = new Line();
            yArrow.Stroke = Brushes.Blue;
            yArrow.StrokeThickness = 4;
            yArrow.StrokeEndLineCap = PenLineCap.Triangle;
            yArrow.X1 = 0;
            yArrow.Y1 = 0;
            yArrow.X2 = 0;
            yArrow.Y2 = 50;
            this.Canvas.Children.Add(yArrow);
        }

        private void DrawShell1dElement(Shell1dElement elem, Brush elemcolor, Brush nodecolor, Double opacity)
        {
            // 要素を描画する
            Polygon polygon = new Polygon();
            polygon.Fill = elemcolor;
            polygon.Stroke = Brushes.Black;
            polygon.StrokeThickness = 1;
            polygon.Opacity = opacity;
            PointCollection Points = new PointCollection();
            for (int i = 0; i < 4; i++)
            {
                Point p = (Point)(elem.Nodes[i].Point * 15);
                p.X += 50;
                p.Y += 500;
                Points.Add(p);

            }
            polygon.Points = Points;
            this.Canvas.Children.Add(polygon);

            // 節点を描画する
            for (int i = 0; i < 4; i++)
            {
                Ellipse ellipse = new Ellipse();
                ellipse.Fill = nodecolor;
                ellipse.Width = 4;
                ellipse.Height = 4;
                Canvas.SetLeft(ellipse, elem.Nodes[i].Point.X * 15 - ellipse.Width * 0.5 + 50);
                Canvas.SetTop(ellipse, elem.Nodes[i].Point.Y * 15 - ellipse.Height * 0.5 + 500);

                this.Canvas.Children.Add(ellipse);
            }
        }
        private void DrawShell2dElement(Shell2dElement elem, Brush elemcolor, Brush nodecolor, Double opacity)
        {
            // 要素を描画する
            Polygon polygon = new Polygon();
            polygon.Fill = elemcolor;
            polygon.Stroke = Brushes.Black;
            polygon.StrokeThickness = 1;
            polygon.Opacity = opacity;
            PointCollection Points = new PointCollection();
            {
                Point p = (Point)(elem.Nodes[0].Point * 15);
                p.X += 50;
                p.Y += 500;
                Points.Add(p);

                p = (Point)(elem.Nodes[4].Point * 15);
                p.X += 50;
                p.Y += 500;
                Points.Add(p);

                p = (Point)(elem.Nodes[1].Point * 15);
                p.X += 50;
                p.Y += 500;
                Points.Add(p);

                p = (Point)(elem.Nodes[5].Point * 15);
                p.X += 50;
                p.Y += 500;
                Points.Add(p);

                p = (Point)(elem.Nodes[2].Point * 15);
                p.X += 50;
                p.Y += 500;
                Points.Add(p);

                p = (Point)(elem.Nodes[6].Point * 15);
                p.X += 50;
                p.Y += 500;
                Points.Add(p);

                p = (Point)(elem.Nodes[3].Point * 15);
                p.X += 50;
                p.Y += 500;
                Points.Add(p);

                p = (Point)(elem.Nodes[7].Point * 15);
                p.X += 50;
                p.Y += 500;
                Points.Add(p);
            }
            polygon.Points = Points;
            this.Canvas.Children.Add(polygon);

            // 節点を描画する
            for (int i = 0; i < 8; i++)
            {
                Ellipse ellipse = new Ellipse();
                ellipse.Fill = nodecolor;
                ellipse.Width = 4;
                ellipse.Height = 4;
                Canvas.SetLeft(ellipse, elem.Nodes[i].Point.X * 15 - ellipse.Width * 0.5 + 50);
                Canvas.SetTop(ellipse, elem.Nodes[i].Point.Y * 15 - ellipse.Height * 0.5 + 500);

                this.Canvas.Children.Add(ellipse);
            }
        }

    }
}

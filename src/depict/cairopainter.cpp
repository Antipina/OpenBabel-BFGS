#include <openbabel/obutil.h>
#include <openbabel/depict/cairopainter.h>

#include <iostream>
using namespace std;

namespace OpenBabel
{

  // Class definition of CairoPainter
  CairoPainter::CairoPainter(int width, int height, string title) : m_surface(0), m_cairo(0),
    m_fontPointSize(12), m_width(width), m_height(height), m_pen_width(1), m_title(title)
  {
  }

  CairoPainter::~CairoPainter()
  {
    if (m_cairo)
      cairo_destroy(m_cairo);
    if (m_surface)
      cairo_surface_destroy(m_surface);
  }

  void CairoPainter::NewCanvas(double width, double height)
  {
    // clean up
    if (m_cairo)
      cairo_destroy(m_cairo);
    if (m_surface)
      cairo_surface_destroy(m_surface);

    // Work out the scaling factor
    double scale_x = m_width / (double) width;
    double scale_y = (m_height-16) / (double) height; // Leave some extra space for the title
    double scale = std::min(scale_x, scale_y);

    // create new surface to paint on
    m_surface = cairo_image_surface_create(CAIRO_FORMAT_ARGB32, static_cast<int> (m_width), static_cast<int> (m_height));
    m_cairo = cairo_create(m_surface);
    cairo_set_source_rgb (m_cairo, 255, 255, 255);
    cairo_paint (m_cairo);
    cairo_set_line_width(m_cairo, m_pen_width);

    // Add the title
    if (!m_title.empty()) {
      this->SetPenColor(OBColor("black"));
      this->SetFontSize(static_cast<int>(16.0));
      OBFontMetrics fm = this->GetFontMetrics(m_title);
      this->DrawText(m_width/2.0 - fm.width/2.0, m_height - fm.height * 0.25, m_title);
    }

    // Translate the over-scaled dimension into the centre
    if (scale < scale_y)
      cairo_translate(m_cairo, 0, m_height/2.0 - scale*height/2.0);
    else
      cairo_translate(m_cairo, m_width/2.0 - scale*width/2.0, 0);
    cairo_scale(m_cairo, scale, scale); // Set a scaling transformation
  }
  
  bool CairoPainter::IsGood() const
  {
    if (!m_cairo)
      return false;
    if (!m_surface)
      return false;
    return true;
  }
      
  void CairoPainter::SetFontSize(int pointSize)
  {
    m_fontPointSize = pointSize;
    cairo_set_font_size(m_cairo, pointSize);
  }

  void CairoPainter::SetFillColor(const OBColor &color)
  {
    cairo_set_source_rgb(m_cairo, color.red, color.green, color.blue);
  }

  void CairoPainter::SetPenColor(const OBColor &color)
  {
    cairo_set_source_rgb(m_cairo, color.red, color.green, color.blue);
  }
      
  void CairoPainter::SetPenWidth(double width)
  {
    m_pen_width = width;
    
  }

  void CairoPainter::DrawLine(double x1, double y1, double x2, double y2)
  {
    cairo_move_to(m_cairo, x1, y1);
    cairo_line_to(m_cairo, x2, y2);
    cairo_stroke(m_cairo);
  }

  void CairoPainter::DrawPolygon(const std::vector<std::pair<double,double> > &points)
  {
    std::vector<std::pair<double,double> >::const_iterator i;
    for (i = points.begin(); i != points.end(); ++i)
      cairo_line_to(m_cairo, i->first, i->second); // note: when called without previous point, 
                                                   //       this function behaves like cairo_move_to 
    cairo_line_to(m_cairo, points.begin()->first, points.begin()->second);
    cairo_fill(m_cairo);
  }

  void CairoPainter::DrawCircle(double x, double y, double r)
  {
    cairo_arc(m_cairo, x, y, r, 0, 2 * M_PI);
    cairo_stroke(m_cairo);
  }

  void CairoPainter::DrawText(double x, double y, const std::string &text)
  {
    cairo_move_to(m_cairo, x, y);
    cairo_show_text(m_cairo, text.c_str());
  }

  OBFontMetrics CairoPainter::GetFontMetrics(const std::string &text)
  {
    cairo_font_extents_t fe;
    cairo_font_extents(m_cairo, &fe);
    cairo_text_extents_t te;
    cairo_text_extents(m_cairo, text.c_str(), &te);

    OBFontMetrics metrics;
    metrics.fontSize = m_fontPointSize;
    metrics.ascent = fe.ascent;
    metrics.descent = -fe.descent;
    metrics.width = te.width;
    metrics.height = fe.height;
    return metrics;
  }
      
  void CairoPainter::WriteImage(const std::string &filename)
  {
    if (!m_cairo || !m_surface)
      return;
    
    cairo_surface_write_to_png(m_surface, filename.c_str());
  }

  static cairo_status_t writeFunction(void* closure, const unsigned char* data, unsigned int length)
  {
    vector<char>* in = reinterpret_cast<vector<char>*>(closure);
    for(int i=0;i<length;++i)
      in->push_back(data[i]);
    return CAIRO_STATUS_SUCCESS;
  }

  void CairoPainter::WriteImage(std::ostream& ofs)
  {
    if (!m_cairo || !m_surface)
      return;
    vector<char> in;
    cairo_surface_write_to_png_stream(m_surface, writeFunction, &in);
    for(int i=0; i<in.size(); ++i)
      ofs << in.at(i);
  }

}


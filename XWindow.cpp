//BL_COPYRIGHT_NOTICE

#include <iostream.h>
#include <stdlib.h>
#include <string.h>

#include "XWindow.H"

static char default_font[] = "8x13";

// draw a character string
void XWindow::putString(int x, int y, const char *str)
{
   int yy = height-y*scale-1;
   XDrawString(display, win, gc, x*scale, yy, str, strlen(str));
   XFlush( display );
}

// clear window
void XWindow::newPage()
{
   XClearWindow( display, win );
   XFlush( display );
}

void XWindow::movePen(int x, int y)
{
   x0 = x*scale;
   y0 = height - y*scale - 1;
}

void XWindow::drawLine(int x, int y)
{

   x *= scale;
   y  = height - y*scale - 1;
   XDrawLine( display, win, gc, x0, y0, x, y );
   XFlush( display );

   x0 = x;
   y0 = y;
}

void XWindow::setLineWidth(int lw) {
  XSetLineAttributes( display, gc, (lw+1)/2,
                      LineSolid, CapRound, JoinRound );
}

void XWindow::flush() {
   XFlush( display );
}

// create a window
XWindow::XWindow(int iwidth, int iheight, const char* str, int iscale)
{
  scale = iscale;
  width = iwidth*scale;
  height = iheight*scale;
  name = new char[strlen(str) + 1];
  strcpy(name,str);
  cmap_colors = NULL;

  screen = 0;
  x0     = 0;
  y0     = 0;


  /* connect to X server */
  if ( (display=XOpenDisplay(NULL)) == NULL ) {
    cerr << "Graphics Display: cannot connect to X server " <<
            XDisplayName(NULL) << "\n";
    abort();
  }

  /* create opaque window with backing store */
  win = XCreateSimpleWindow(display,
	                    RootWindow(display,screen),
			    100, 100, width, height, 4,
			    BlackPixel(display,screen),
			    WhitePixel(display,screen) );

  /* initialize size hint property for window manager */
  XSizeHints size_hints;
  size_hints.flags = PSize | PMaxSize | PMinSize;
  size_hints.width = (width+1);
  size_hints.height = (height+1);
  size_hints.max_width = size_hints.min_width = (width+1);
  size_hints.max_height = size_hints.min_height = (height+1);

  /* set properties for window manager (always before mapping) */
  XSetStandardProperties( display, win, name,
			  name,
			  None, NULL, 0, &size_hints );

  /* set window manager hints */
  XWMHints	wm_hints;
  wm_hints.initial_state = NormalState;
  wm_hints.input = False;
  wm_hints.flags = StateHint | InputHint;

  XSetWMHints(display, win, &wm_hints);

  /* set the class hints */
  XClassHint class_hints;
  class_hints.res_name = name;
  class_hints.res_class = name;

  XSetClassHint(display, win, &class_hints);

  /* specify default font */
  if ((the_font = XLoadQueryFont(display,default_font)) == NULL) {
     cerr << "Cannot open font: " << default_font << "\n";
     abort();
  };

  /* create default graphics context */
  XGCValues gcvalues;
  gcvalues.font = the_font->fid;

  gc = XCreateGC( display, win, GCFont, &gcvalues );

  /* specify black foreground since default may be white on white */
  XSetForeground( display, gc, BlackPixel(display,screen) );

  /* set line attributes */
  XSetLineAttributes( display, gc, 1, LineSolid,
                      CapRound, JoinRound );

  /* create a colormap */
  nplanes = DisplayPlanes(display,screen);
  if (nplanes > 2) {
     cmap = DefaultColormap(display,screen);
     XSetWindowColormap(display,win,cmap);
  }; /* psuedocolor display */

  /* select event types wanted */
  // XSelectInput( display, win, 0L );
  XSelectInput( display, win, ExposureMask | StructureNotifyMask );
  /* display window */
  XMapWindow( display, win );
  // start  an  event loop that runs  until the window has received its
  // initial exposure events.
  while(1) {
    XEvent report;
    XNextEvent(display, &report);
    switch(report.type) {
    case Expose:
	if (report.xexpose.count != 0)
	    break;
	goto op_end;
    default:
	break;
    }
     }
  op_end:XSelectInput( display, win, 0L );
  /* some displays allow backing store,  others don't */
  if (DoesBackingStore(DefaultScreenOfDisplay(display))) {
      XSetWindowAttributes setwinattr;
      setwinattr.backing_store = Always;
      XChangeWindowAttributes(display, win,
			      CWBackingStore, &setwinattr);
  }

}

// remove a window
XWindow::~XWindow()
{
   // how to free colors?
  XFreeFont(display, the_font);
   XFreeGC( display, gc );
   XCloseDisplay( display );
   delete name;
   name = NULL;
   delete cmap_colors;
   width = height = -1;
}

// set font
void  XWindow::setFont(char *font)
{
  XGCValues gcvalues;
  if ((the_font = XLoadQueryFont(display,font)) == NULL) {
     cerr << "Cannot open font: " << font << "\n";
     abort();
  };
  gcvalues.font = the_font->fid;
  XChangeGC( display, gc, GCFont, &gcvalues );
}


void XWindow::drawBox(int x1, int y1, int x2, int y2)
{
   x1 *= scale;
   x2 *= scale;
   y1  = height - y1*scale - 1;
   y2  = height - y2*scale - 1;
   XDrawRectangle(display,win,gc,x1,y1,x2-x1,y2-y1);
   XFlush( display );

   x0 = x1;
   y0 = y1;
}

void XWindow::setfgColor(const char *colorname) {

  // if monotone display, NoOp.
  if (nplanes <= 2) return;

  XColor color, exact;
  if (!XAllocNamedColor(display,cmap,colorname,&color,&exact)) {
     cerr << "cant alloc colortable cell for color = " << colorname << endl;
  };
  color.flags = DoRed | DoGreen | DoBlue;
  XSetForeground( display, gc, color.pixel );
}

void XWindow::setbgColor(const char *colorname) {

  // if monotone display, NoOp.
  if (nplanes <= 2) return;

  XColor color, exact;
  if (!XAllocNamedColor(display,cmap,colorname,&color,&exact)) {
     cerr << "cant alloc colortable cell for color = " << colorname << endl;
  };
  XSetWindowBackground( display, win, color.pixel );
  XClearWindow(display, win);
}

void XWindow::setfgColor(int col) {
  // if monotone display, NoOp.
  if (nplanes <= 2) return;

  if ( (col < 0) || (col >= num_cmap_colors) ) {
     cout << "color out of range: " << col << endl;
     return;
  };
  XSetForeground( display, gc, cmap_colors[col]);
}

void XWindow::setbgColor(int col) {
  // if monotone display, NoOp.
  if (nplanes <= 2) return;

  if ( (col < 0) || (col >= num_cmap_colors) ) {
     cout << "color out of range: " << col << endl;
     return;
  };
  XSetWindowBackground( display, win, cmap_colors[col]);
  XClearWindow(display, win);
}

void XWindow::defineCmap(unsigned short *red, unsigned short *green,
                         unsigned short *blue, int num) {

  // if monotone display, do nothing
  if (nplanes <= 2)  return;

  XColor *cols = new XColor[num];
  cmap_colors = new unsigned long[num];
  if (!XAllocColorCells(display,cmap,True,NULL,0,cmap_colors,num)) {
     cerr << "could not alloc " << num <<" colors in colortable" << endl;
     abort();
  };
  num_cmap_colors = num;
  for (int i = 0; i < num; i++) {
     cols[i].red = red[i];
     cols[i].blue = blue[i];
     cols[i].green = green[i];
     cols[i].flags = DoRed | DoGreen | DoBlue;
     cols[i].pixel = cmap_colors[i];
  };
  XStoreColors(display,cmap,cols,num);
  delete cols;
}

int XWindow::getMouse(int &x, int &y)
{
  XEvent report;
  int    ib;

  /* flush event queue */
  while ( XCheckMaskEvent( display, ~0L, &report ) );

  /* select button press events */
  XSelectInput( display, win, ButtonPressMask );

  /* loop until a button press event is detected */
  ib = 0;
  while (ib == 0) {
    XNextEvent( display, &report );
    if (report.type == ButtonPress) {
      ib = report.xbutton.button;
      x = report.xbutton.x;
      y = report.xbutton.y;
    };
  };

  x /= scale;
  y  = height - y - 1;
  y /= scale;

  /* deselect button press events */
  XSelectInput( display, win, 0L );

  return ib;
}



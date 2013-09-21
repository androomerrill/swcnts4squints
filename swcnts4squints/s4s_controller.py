from PyQt4.QtGui import QApplication

from s4s_view import S4SView

class S4SController(object):

    def __init__(self, args, model=None):

        app = QApplication(args)
        self.model = model

        self.view = S4SView(self, self.model)
        self.model.init()
        self.view.show()
        app.exec_()

    def refresh_view(self):
        self.view.update_app_view()

from flask_wtf import Form
from wtforms import StringField, BooleanField, TextAreaField
from wtforms.validators import DataRequired

class LoginForm(Form):
        openid = StringField('openid', validators=[DataRequired()])
        remember_me = BooleanField('remember_me', default=False)

class ProteinInputForm(Form):
        accessionInput = TextAreaField(u'Accession Values', validators=[DataRequired()])

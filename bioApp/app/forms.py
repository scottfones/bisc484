from flask_wtf import Form
from flask_wtf.file import FileField, FileRequired, FileAllowed
from wtforms.validators import DataRequired
from wtforms import StringField, BooleanField, TextAreaField, validators

class LoginForm(Form):
    openid = StringField('openid', validators=[DataRequired()])
    remember_me = BooleanField('remember_me', default=False)

class ProteinInputForm(Form):
    accessionInput = TextAreaField(u'Accession Values', validators=[DataRequired()])

class ABIInputForm(Form):
    abi_file =  FileField(u'ABI File', validators=
                          [FileAllowed(['abi', 'ab1', 'ABI', 'AB1'], 'ABI Files Only')])

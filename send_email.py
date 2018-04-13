# We start by only importing only the classes we need, this also saves us from
# having to use the full module name later.
from email.MIMEMultipart import MIMEMultipart
from email.MIMEText import MIMEText
# Then we compose some of the basic message headers:
fromaddr = "sender@gmail.com"
toaddr = "receiver@gmail.com"
msg = MIMEMultipart()
msg['From'] = fromaddr
msg['To'] = toaddr
msg['Subject'] = "Python email"
# Next, we attach the body of the email to the MIME message:
body = "Python test mail"
msg.attach(MIMEText(body, 'plain'))
# For sending the mail, we have to convert the object to a string, and then
# use the same prodecure as above to send using the SMTP server..
import smtplib
server = smtplib.SMTP('smtp.gmail.com', 587)
server.ehlo()
server.starttls()
server.ehlo()
server.login("some_email", "PASSWORD")
text = msg.as_string()
server.sendmail(fromaddr, toaddr, text)

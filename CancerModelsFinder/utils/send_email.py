import smtplib
import os


def send(body, subject):
   # TODO changed to local
   sender_email = 'no-reply@ebi.ac.uk'
   port = 25
   smtp_server = "smtp.ebi.ac.uk"
   receiver_emails = ['tushar@ebi.ac.uk']

   email_text = """\
From: %s
To: %s
Subject: %s

%s
""" % (sender_email, receiver_emails, subject, body)

   try:
       server = smtplib.SMTP(smtp_server, 25)
       server.connect(smtp_server, port)
       server.ehlo()
       server.sendmail(sender_email, receiver_emails, email_text)
       server.quit()
       print("Email sent successfully!")
   except Exception as ex:
       print("Something went wrong while sending mail….", ex)


#send('test body', 'test subject for tushar')
